!Fitting Hd is a nonlinear least square fit problem, 2-step method is my favourite:
!    Step 1: A fast but inexact hopping to explore the phase space
!        Here we adopt the iterative pseudolinear equations specially designed for SurfGen
!        At a given geometry, diagonalize the current Hd for adiabatic states {phi_i}
!        Fix {phi_i}, solve the linear least square fit equations to update Hd
!    Step 2: A rigorous local minimizer based on the best estimation explored
!
!Nomenclature:
!The Lagrangian is defined as:
!    L = sum( rho^2 || H_ad^d - H_ad^ab ||_F^2 + || ▽H_ad^d - ▽H_ad^ab ||_F^2, over point)
!      + sum( rho^2 || H_nd^d - H_nd^ab ||_F^2 + || ▽H_nd^d - ▽H_nd^ab ||_F^2, over DegeneratePoint)
!      + sum( rho^2 || H_ad^d - H_ad^ab ||_F^2, over ArtifactPoint) + 
!where rho is the unit converter from energy to energy gradient (HdLSF_EnergyScale),
!    subscript ad means adiabatic representation, nd means nondegenerate representation, F means Frobenius norm
!    superscript d means diabatz, ab means ab initio
!    definition of point, DegeneratePoint, ArtifactPoint see module Basic
!We may also add a regularization to Lagrangian, namely tikhonov regularization: tau * || c ||_2^2
!where tau is the parameterized KKT multiplier (HdLSF_Regularization), c is undetermined parameter vector in Hd
!
!Implementation detail:
!Here we sort y by LSFPoint -> ArtifactPoint, each point provides Hamiltonian column by column first,
!    then gradients by istate -> jstate ( >= istate ) from 1st direction to the last
!c is sorted by istate -> jstate ( >= istate ) -> iorder according to module DiabaticHamiltonian
module HdLeastSquareFit
    use Basic
    use DiabaticHamiltonian
    implicit none

!Parameter
    !General solver control:
        real*8::HdLSF_Regularization=0d0!Instead of solving the KKT multiplier, let it be a parameter
        !Available solvers: pseudolinear, TrustRegion, LineSearch
        !Choose the nonlinear optimization solver: a single solver or a 2-step solver
        !    pseudolinear_X combining pseudolinear and another solver X
        !If you see insufficient memory, use only LineSearch (and LineSearcher = LBFGS or ConjugateGradient)
        !Warning: pseudolinear is a quick hopper by vanishing the linear part of the gradient
        !    However, it is not necessarily able to solve the fitting alone, 
        !    unless the minimum coincidentally also has the nonlinear part of the gradient = 0
        character*32::HdLSF_Solver='pseudolinear_TrustRegion'
        !Max ineration control. Hopper = pseudolinear. LocalMinimizer = TrustRegion, LineSearch
        integer::HdLSF_MaxHopperIteration=100,HdLSF_MaxLocalMinimizerIteration=1000,HdLSF_Max2StepIteration=10
    !pseudolinear:
        integer::HdLSF_pseudolinearFollowFreq=1,&!Every how many steps print fitting progress
            HdLSF_pseudolinearMaxMonotonicalIncrease=10!Terminate after how many monotonically increasing iterations
    !LineSearch:
        character*32::HdLSF_LineSearcher='ConjugateGradient'!Available: LBFGS, ConjugateGradient
        logical::UseStrongWolfe=.true.!Whether use strong Wolfe condition instead of Wolfe condition
        !LBFGS:
            integer::HdLSF_LBFGSMemory=10!Memory usage control, [3,30] is recommended (must > 0)
        !ConjugateGradient:
            character*2::HdLSF_ConjugateGradientSolver='DY'!Available: DY (Dai-Yun), PR (Polak-Ribiere+)

!HdLeastSquareFit module only variable
    integer::HdLSF_NData!Number of fitting data
    real*8::HdLSF_EnergyScale!Scale energy residue in Lagrangian for its unit difference from gradient
    real*8::HdLSF_SqrtRegularization,HdLSF_EnergyScaleSquare
    !Work space
        real*8,allocatable,dimension(:)::HdLSF_energy
        !dc = ▽_c. H, dH, dcH, dcdH are in representation of the data point 
        !For dcH & dcdH, ▽_c only operates on H & ▽H
        !phi is the basis matrix of the data point representation in diabatic representation
        real*8,allocatable,dimension(:,:)::HdLSF_H,HdLSF_phi
        real*8,allocatable,dimension(:,:,:)::HdLSF_dH,HdLSF_dcH
        real*8,allocatable,dimension(:,:,:,:)::HdLSF_dcdH
        real*8,allocatable,dimension(:)::HdLSF_f!Stores expansion basis function values
        real*8,allocatable,dimension(:,:)::HdLSF_fd!Stores expansion basis function gradient values
        !Pseudolinear: W is diagonal so only store its diagonal vector, MT = M^T, NT = N^T
            real*8,allocatable,dimension(:)::HdLSF_y,HdLSF_W
            real*8,allocatable,dimension(:,:)::HdLSF_M,HdLSF_MT
        !Local minimizer
            !dHd is in diabatic representaton
            !For description of dcphi, see Nonadiabatic.deigvec_ByKnowneigval_dA (with ▽ replaced by ▽_c)
            !For dcHrep & dcdHrep, ▽_c also operates on the representation basis
            real*8,allocatable,dimension(:,:,:)::HdLSF_dHd,HdLSF_dcphi,HdLSF_dcHrep
            real*8,allocatable,dimension(:,:,:,:)::HdLSF_dcdHrep
            !Trust region
                real*8,allocatable,dimension(:)::HdLSF_spResidue,HdLSF_sdegpResidue!Residue of a single (degenerate) data point
                real*8,allocatable,dimension(:,:)::HdLSF_spJacobian,HdLSF_sdegpJacobian!Jacobian of a single (degenerate) data point

contains
!The initializer for HdLeastSquareFit module
subroutine InitializeHdLeastSquareFit()
    integer::ip,istate,jstate
    real*8::MaxEnergy,MaxGrad,temp
    !A data point provides NStates adiabatic energies, InternalDimension x NStates x NStates ▽H
    DataPerPoint=NStates+NStates*(NStates+1)/2*InternalDimension!Upper triangle is redundant
    !A degenerate data point provides NStates order H, InternalDimension x NStates x NStates ▽H
    DataPerDegeneratePoint=NStates*(NStates+1)/2*(InternalDimension+1)!Upper triangle is redundant
    MaxEnergy=0d0
    MaxGrad=0d0
    do ip=1,NPoints
        temp=maxval(Abs(point(ip).energy))
        if(temp>MaxEnergy) MaxEnergy=temp
        do istate=1,NStates
            do jstate=istate,NStates
                temp=maxval(abs(point(ip).dH(:,jstate,istate)))
                if(temp>MaxGrad) MaxGrad=temp
            end do
        end do
    end do
    HdLSF_EnergyScale=MaxGrad/MaxEnergy
    HdLSF_EnergyScaleSquare=HdLSF_EnergyScale*HdLSF_EnergyScale
    HdLSF_SqrtRegularization=dSqrt(HdLSF_Regularization)
end subroutine InitializeHdLeastSquareFit

subroutine FitHd()!Fit Hd with the designated solver
    integer::istate,jstate
    real*8,allocatable,dimension(:)::c
    !Initialize
        HdLSF_NData=DataPerPoint*NPoints+DataPerDegeneratePoint*NDegeneratePoints+NStates*NArtifactPoints
        !Initial value of c
        allocate(c(NExpansionCoefficients))
        call HdEC2c(HdEC,c,NExpansionCoefficients)
        !Allocate global work space
            if(allocated(HdLSF_energy)) deallocate(HdLSF_energy)
            allocate(HdLSF_energy(NStates))
            if(allocated(HdLSF_H)) deallocate(HdLSF_H)
            allocate(HdLSF_H(NStates,NStates))
            if(allocated(HdLSF_phi)) deallocate(HdLSF_phi)
            allocate(HdLSF_phi(NStates,NStates))
            if(allocated(HdLSF_dH)) deallocate(HdLSF_dH)
            allocate(HdLSF_dH(InternalDimension,NStates,NStates))
            if(allocated(HdLSF_dcH)) deallocate(HdLSF_dcH)
            allocate(HdLSF_dcH(NExpansionCoefficients,NStates,NStates))
            if(allocated(HdLSF_dcdH)) deallocate(HdLSF_dcdH)
            allocate(HdLSF_dcdH(NExpansionCoefficients,InternalDimension,NStates,NStates))
            if(allocated(HdLSF_f)) deallocate(HdLSF_f)
            allocate(HdLSF_f(NExpansionBasis))
            if(allocated(HdLSF_fd)) deallocate(HdLSF_fd)
            allocate(HdLSF_fd(InternalDimension,NExpansionBasis))
    !Fit Hd
    select case(HdLSF_Solver)
        !Single solvers
        case('pseudolinear')
            call pseudolinear(c)
        case('TrustRegion')
            call TrustRegionInterface(c)
        case('LBFGS')
            call LBFGSInterface(c)
        case('ConjugateGradient')
            call ConjugateGradientInterface(c)
        !2-step solvers
        case('pseudolinear_TrustRegion')
            do istate=1,HdLSF_Max2StepIteration
                call pseudolinear(c)
                call TrustRegionInterface(c)
            end do
        case('pseudolinear_LBFGS')
            do istate=1,HdLSF_Max2StepIteration
                call pseudolinear(c)
                call LBFGSInterface(c)
            end do
        case('pseudolinear_ConjugateGradient')
            do istate=1,HdLSF_Max2StepIteration
                call pseudolinear(c)
                call ConjugateGradientInterface(c)
            end do
        case default!Throw a warning
            write(*,'(1x,A50,1x,A32)')'Program abort: unsupported least square fit solver',HdLSF_Solver
            stop
    end select
    !Clean up
        deallocate(c)
        !Global work space
            deallocate(HdLSF_energy)
            deallocate(HdLSF_H)
            deallocate(HdLSF_phi)
            deallocate(HdLSF_dH)
            deallocate(HdLSF_dcH)
            deallocate(HdLSF_dcdH)
            deallocate(HdLSF_f)
            deallocate(HdLSF_fd)
end subroutine FitHd

!--------------- Solvers ---------------
    subroutine pseudolinear(cmin)!Fit Hd by pseudolinear method
        real*8,dimension(NExpansionCoefficients),intent(inout)::cmin
        integer::indice,ip,i
        real*8::cchange,L,Lold,Lmin,RMSDenergy,RMSDdH,RMSDDegH,RMSDDegdH
        real*8,allocatable,dimension(:)::b,c
        real*8,allocatable,dimension(:,:)::A
        !Initialize
            !Allocate global pseudolinear work space
                if(allocated(HdLSF_y)) deallocate(HdLSF_y)
                allocate(HdLSF_y(HdLSF_NData))
                if(allocated(HdLSF_W)) deallocate(HdLSF_W)
                allocate(HdLSF_W(HdLSF_NData))
                !W will not change through out the solving procedure, so fill in its value now
                    indice=1
                    do ip=1,NPoints
                        HdLSF_W(indice:indice+NStates-1)=point(ip).weight!Energy
                        indice=indice+NStates
                        do i=1,NStates!▽H
                            HdLSF_W(indice:indice+InternalDimension-1)=point(ip).weight!Diagonal
                            indice=indice+InternalDimension
                            HdLSF_W(indice:indice+InternalDimension*(NStates-i)-1)=2d0*point(ip).weight!Off-diagonal
                            indice=indice+InternalDimension*(NStates-i)
                        end do
                    end do
                    do ip=1,NDegeneratePoints
                        do i=1,NStates!H
                            HdLSF_W(indice)=DegeneratePoint(ip).weight!Diagonal
                            HdLSF_W(indice+1:indice+NStates-i)=2d0*DegeneratePoint(ip).weight!Off-diagonal
                            indice=indice+NStates-i+1
                        end do
                        do i=1,NStates!▽H
                            HdLSF_W(indice:indice+InternalDimension-1)=DegeneratePoint(ip).weight!Diagonal
                            indice=indice+InternalDimension
                            HdLSF_W(indice:indice+InternalDimension*(NStates-i)-1)=2d0*DegeneratePoint(ip).weight!Off-diagonal
                            indice=indice+InternalDimension*(NStates-i)
                        end do
                    end do
                    do ip=1,NArtifactPoints
                        HdLSF_W(indice:indice+NStates-1)=ArtifactPoint(ip).weight!Energy only
                        indice=indice+NStates
                    end do
                if(allocated(HdLSF_M)) deallocate(HdLSF_M)
                allocate(HdLSF_M(NExpansionCoefficients,HdLSF_NData))
                if(allocated(HdLSF_MT)) deallocate(HdLSF_MT)
                allocate(HdLSF_MT(HdLSF_NData,NExpansionCoefficients))
            !Initialize linear least square fit and Lagrangian minimum
                allocate(c(NExpansionCoefficients))
                c=cmin
                allocate(b(NExpansionCoefficients))
                b=c
                allocate(A(NExpansionCoefficients,NExpansionCoefficients))
                call LSFMatrices_L(A,b,Lmin)
                L=Lmin
                indice=0
        !Solve
        call showtime()
        write(*,'(1x,A43)')'Explore phase space by pseudolinear hopping'
        do i=1,HdLSF_MaxHopperIteration!Main loop
            Lold=L!Prepare
            call My_dposv(A,b,NExpansionCoefficients)
            cchange=dot_product(b-c,b-c)
            if(cchange<1d-30) then!Convergence standard: || c_new - c_old ||^2 < 1d-30
                call L_RMSD(b,L,RMSDenergy,RMSDdH,RMSDDegH,RMSDDegdH)
                if(L<Lmin) then
                    Lmin=L
                    cmin=b
                end if
                write(*,*)
                write(*,*)'Hd expansion coefficients have converged at iteration',i
                write(*,*)'Lagrangian =',L
                write(*,*)'Lowest Lagrangian encountered =',Lmin
                write(*,*)'RMSD over regular data points:'
                write(*,*)'     E =',RMSDenergy/cm_1InAu,'cm-1'
                write(*,*)'    dH =',RMSDdH,'a.u.'
                if(NDegeneratePoints>0) then
                    write(*,*)'RMSD over almost degenerate data points:'
                    write(*,*)'     H =',RMSDDegH/cm_1InAu,'cm-1'
                    write(*,*)'    dH =',RMSDDegdH,'a.u.'
                end if
                exit
            else if(mod(i,HdLSF_pseudolinearFollowFreq)==0) then
                c=b
                call LSFMatrices_L_RMSD(A,b,L,RMSDenergy,RMSDdH,RMSDDegH,RMSDDegdH)
                if(L<Lmin) then
                    Lmin=L
                    cmin=c
                end if
                if(L>Lold) then
                    if(indice>HdLSF_pseudolinearMaxMonotonicalIncrease) then
                        write(*,'(1x,A117)')'Pseudolinear hopper warning: hopping is not making progress, but Hd expansion coefficients have not met accuracy goal'
                        exit
                    else
                        indice=indice+1
                    end if
                else
                    indice=0
                end if
                call ShowTime()
                write(*,*)'Iteration',i
                write(*,*)'Change of Hd expansion coefficients =',cchange
                write(*,*)'Lagrangian =',L
                write(*,*)'Lowest Lagrangian encountered =',Lmin
                write(*,*)'RMSD over regular data points:'
                write(*,*)'     E =',RMSDenergy/cm_1InAu,'cm-1'
                write(*,*)'    dH =',RMSDdH,'a.u.'
                if(NDegeneratePoints>0) then
                    write(*,*)'RMSD over almost degenerate data points:'
                    write(*,*)'     H =',RMSDDegH/cm_1InAu,'cm-1'
                    write(*,*)'    dH =',RMSDDegdH,'a.u.'
                end if
            else
                c=b
                call LSFMatrices_L(A,b,L)
                if(L<Lmin) then
                    Lmin=L
                    cmin=c
                end if
                if(L>Lold) then
                    if(indice>HdLSF_pseudolinearMaxMonotonicalIncrease) then
                        write(*,'(1x,A117)')'Pseudolinear hopper warning: hopping is not making progress, but Hd expansion coefficients have not met accuracy goal'
                        exit
                    else
                        indice=indice+1
                    end if
                else
                    indice=0
                end if
            end if
        end do
        !Clean up
            !Local work space
                deallocate(b)
                deallocate(c)
                deallocate(A)
            !Global pseudolinear work space
                deallocate(HdLSF_y)
                deallocate(HdLSF_W)
                deallocate(HdLSF_M)
                deallocate(HdLSF_MT)
        !Output
        if(i>HdLSF_MaxHopperIteration) then!Throw a warning
            write(*,*)
            write(*,'(1x,A52)')'Failed pseudolinear hopping: max iteration exceeded!'
        end if
        call L_RMSD(cmin,L,RMSDenergy,RMSDdH,RMSDDegH,RMSDDegdH)
        write(*,'(1x,A40)')'Best estimation of pseudolinear hopping:'
        write(*,*)'Lagrangian =',L
        write(*,*)'RMSD over regular data points:'
        write(*,*)'     E =',RMSDenergy/cm_1InAu,'cm-1'
        write(*,*)'    dH =',RMSDdH,'a.u.'
        if(NDegeneratePoints>0) then
            write(*,*)'RMSD over almost degenerate data points:'
            write(*,*)'     H =',RMSDDegH/cm_1InAu,'cm-1'
            write(*,*)'    dH =',RMSDDegdH,'a.u.'
        end if
        write(*,'(1x,A30)')'Save Hd expansion coefficients'
        call c2HdEC(cmin,HdEC,NExpansionCoefficients)
        call WriteHdExpansionCoefficients(HdEC,NStates,NOrder)
        contains
            !The general form of weighted linear least square fit with pseudoregularization is:
            !    A c = b, where A = M . W . M^T + tau, b = M . W . y, c & tau have been explained at header,
            !    y is the data vector, M^T . c is the fitting prediction of y, W is the weight
            !To save memory, off-diagonals are treated as twice weighed
            !Input:  b = current c
            !Output: A harvests A, b harvests b, L harvests Lagrangian
            subroutine LSFMatrices_L(A,b,L)
                real*8,dimension(NExpansionCoefficients,NExpansionCoefficients),intent(out)::A
                real*8,dimension(NExpansionCoefficients),intent(inout)::b
                real*8,intent(out)::L
                integer::ip,istate,jstate,indicerow
                real*8::Ltemp
                !Initialize
                    L=HdLSF_Regularization*dot_product(b,b)!Regularization
                    call c2HdEC(b,HdEC,NExpansionCoefficients)
                !Construct M^T and y, add least square fit penalty to Lagrangian
                indicerow=1!Start from 1st row
                do ip=1,NPoints!Regular data points
                    call AdiabaticEnergy_dH_State_f_fd(point(ip).geom,HdLSF_energy,HdLSF_dH,HdLSF_phi,HdLSF_f,HdLSF_fd)!Adiabatic representation
                    call dAssignBasisPhaseBydH(HdLSF_phi,HdLSF_dH,point(ip).dH,Ltemp,InternalDimension,NStates)!Assign basis phase
                    !Energy
                    HdLSF_energy=HdLSF_energy-point(ip).energy
                    L=L+point(ip).weight*(Ltemp+HdLSF_EnergyScaleSquare*dot_product(HdLSF_energy,HdLSF_energy))
                    HdLSF_dcH=sy3UnitaryTransformation(dcHd_ByKnownf(HdLSF_f),HdLSF_phi,NExpansionCoefficients,NStates)
                    HdLSF_y(indicerow:indicerow+NStates-1)=HdLSF_EnergyScaleSquare*point(ip).energy
                    forall(istate=1:NStates)
                        HdLSF_MT(indicerow+istate-1,:)=HdLSF_EnergyScaleSquare*HdLSF_dcH(:,istate,istate)
                    end forall
                    indicerow=indicerow+NStates
                    !▽H (▽_c phi is neglected)
                    HdLSF_dcdH=sy4UnitaryTransformation(dcdHd_ByKnownfdT(transpose(HdLSF_fd)),HdLSF_phi,NExpansionCoefficients,InternalDimension,NStates)
                    do istate=1,NStates
                        do jstate=istate,NStates
                            HdLSF_y( indicerow:indicerow+InternalDimension-1)=point(ip).dH(:,jstate,istate)
                            HdLSF_MT(indicerow:indicerow+InternalDimension-1,:)=transpose(HdLSF_dcdH(:,:,jstate,istate))
                            indicerow=indicerow+InternalDimension
                        end do
                    end do
                end do
                do ip=1,NDegeneratePoints!Almost degenerate data points
                    call NondegenerateH_dH_State_f_fd(DegeneratePoint(ip).geom,HdLSF_H,HdLSF_dH,HdLSF_phi,HdLSF_f,HdLSF_fd)!Nondegenerate representation
                    call dFixHPhase_AssignBasisPhaseBydH(HdLSF_H,HdLSF_phi,HdLSF_dH,DegeneratePoint(ip).dH,Ltemp,InternalDimension,NStates)!Assign basis phase
                    !H
                    forall(istate=1:NStates,jstate=1:NStates,istate>=jstate)
                        HdLSF_H(istate,jstate)=HdLSF_H(istate,jstate)-DegeneratePoint(ip).H(istate,jstate)
                    end forall
                    L=L+DegeneratePoint(ip).weight*(Ltemp+HdLSF_EnergyScaleSquare*dsyFrobeniusSquare(HdLSF_H,NStates))
                    HdLSF_dcH=sy3UnitaryTransformation(dcHd_ByKnownf(HdLSF_f),HdLSF_phi,NExpansionCoefficients,NStates)
                    do istate=1,NStates
                        HdLSF_y( indicerow:indicerow+NStates-istate)=HdLSF_EnergyScaleSquare*DegeneratePoint(ip).H(istate:NStates,istate)
                        HdLSF_MT(indicerow:indicerow+NStates-istate,:)=HdLSF_EnergyScaleSquare*transpose(HdLSF_dcH(:,istate:NStates,istate))
                        indicerow=indicerow+NStates-istate+1
                    end do
                    !▽H (▽_c phi is neglected)
                    HdLSF_dcdH=sy4UnitaryTransformation(dcdHd_ByKnownfdT(transpose(HdLSF_fd)),HdLSF_phi,NExpansionCoefficients,InternalDimension,NStates)
                    do istate=1,NStates
                        do jstate=istate,NStates
                            HdLSF_y( indicerow:indicerow+InternalDimension-1)=DegeneratePoint(ip).dH(:,jstate,istate)
                            HdLSF_MT(indicerow:indicerow+InternalDimension-1,:)=transpose(HdLSF_dcdH(:,:,jstate,istate))
                            indicerow=indicerow+InternalDimension
                        end do
                    end do
                end do
                Ltemp=0d0
                do ip=1,NArtifactPoints!Unreliable data points, energy only
                    call AdiabaticEnergy_State_f(ArtifactPoint(ip).geom,HdLSF_energy,HdLSF_phi,HdLSF_f)
                    HdLSF_energy=HdLSF_energy-ArtifactPoint(ip).energy
                    Ltemp=Ltemp+ArtifactPoint(ip).weight*dot_product(HdLSF_energy,HdLSF_energy)
                    HdLSF_dcH=sy3UnitaryTransformation(dcHd_ByKnownf(HdLSF_f),HdLSF_phi,NExpansionCoefficients,NStates)
                    HdLSF_y(indicerow:indicerow+NStates-1)=HdLSF_EnergyScaleSquare*ArtifactPoint(ip).energy
                    forall(istate=1:NStates)
                        HdLSF_MT(indicerow+istate-1,:)=HdLSF_EnergyScaleSquare*HdLSF_dcH(:,istate,istate)
                    end forall
                    indicerow=indicerow+NStates
                end do
                L=L+HdLSF_EnergyScaleSquare*Ltemp
                !Done construction, put them into A and b
                HdLSF_M=transpose(HdLSF_MT)
                forall(ip=1:HdLSF_NData)
                    HdLSF_M(:,ip)=HdLSF_M(:,ip)*HdLSF_W(ip)
                end forall
                A=matmul(HdLSF_M,HdLSF_MT)
                forall(ip=1:NExpansionCoefficients)!Regularization
                    A(ip,ip)=A(ip,ip)+HdLSF_Regularization
                end forall
                b=matmul(HdLSF_M,HdLSF_y)
            end subroutine LSFMatrices_L
            !Input:  b = current c
            !Output: A harvests A, b harvests b, L harvests Lagrangian
            !        RMSDenergy/dH harvests root mean square deviation of adiabatic energy/dH over point,
            !        RMSDDegH/dH harvests root mean square deviation of nondegenerate H/dH over DegeneratePoint
            subroutine LSFMatrices_L_RMSD(A,b,L,RMSDenergy,RMSDdH,RMSDDegH,RMSDDegdH)
                real*8,dimension(NExpansionCoefficients,NExpansionCoefficients),intent(out)::A
                real*8,dimension(NExpansionCoefficients),intent(inout)::b
                real*8,intent(out)::L,RMSDenergy,RMSDdH,RMSDDegH,RMSDDegdH
                integer::ip,istate,jstate,indicerow
                real*8::Ltemp,temp
                !Initialize
                    L=HdLSF_Regularization*dot_product(b,b)!Regularization
                    call c2HdEC(b,HdEC,NExpansionCoefficients)
                !Construct M^T and y, add least square fit penalty to Lagrangian
                indicerow=1!Start from 1st row
                RMSDenergy=0d0
                RMSDdH=0d0
                do ip=1,NPoints!Regular data points, compute RMSD
                    call AdiabaticEnergy_dH_State_f_fd(point(ip).geom,HdLSF_energy,HdLSF_dH,HdLSF_phi,HdLSF_f,HdLSF_fd)!Adiabatic representation
                    call dAssignBasisPhaseBydH(HdLSF_phi,HdLSF_dH,point(ip).dH,Ltemp,InternalDimension,NStates)!Assign basis phase
                    RMSDdH=RMSDdH+Ltemp
                    !Energy
                    HdLSF_energy=HdLSF_energy-point(ip).energy
                    temp=dot_product(HdLSF_energy,HdLSF_energy)
                    RMSDenergy=RMSDenergy+temp
                    L=L+point(ip).weight*(Ltemp+HdLSF_EnergyScaleSquare*temp)
                    HdLSF_dcH=sy3UnitaryTransformation(dcHd_ByKnownf(HdLSF_f),HdLSF_phi,NExpansionCoefficients,NStates)
                    HdLSF_y(indicerow:indicerow+NStates-1)=HdLSF_EnergyScaleSquare*point(ip).energy
                    forall(istate=1:NStates)
                        HdLSF_MT(indicerow+istate-1,:)=HdLSF_EnergyScaleSquare*HdLSF_dcH(:,istate,istate)
                    end forall
                    indicerow=indicerow+NStates
                    !▽H (▽_c phi is neglected)
                    HdLSF_dcdH=sy4UnitaryTransformation(dcdHd_ByKnownfdT(transpose(HdLSF_fd)),HdLSF_phi,NExpansionCoefficients,InternalDimension,NStates)
                    do istate=1,NStates
                        do jstate=istate,NStates
                            HdLSF_y( indicerow:indicerow+InternalDimension-1)=point(ip).dH(:,jstate,istate)
                            HdLSF_MT(indicerow:indicerow+InternalDimension-1,:)=transpose(HdLSF_dcdH(:,:,jstate,istate))
                            indicerow=indicerow+InternalDimension
                        end do
                    end do
                end do
                RMSDenergy=dSqrt(RMSDenergy/NStates/NPoints)
                RMSDdH=dSqrt(RMSDdH/(InternalDimension*NStates*NStates)/NPoints)
                RMSDDegH=0d0
                RMSDDegdH=0d0
                do ip=1,NDegeneratePoints!Almost degenerate data points, compute RMSDDeg
                    call NondegenerateH_dH_State_f_fd(DegeneratePoint(ip).geom,HdLSF_H,HdLSF_dH,HdLSF_phi,HdLSF_f,HdLSF_fd)!Nondegenerate representation
                    call dFixHPhase_AssignBasisPhaseBydH(HdLSF_H,HdLSF_phi,HdLSF_dH,DegeneratePoint(ip).dH,Ltemp,InternalDimension,NStates)!Assign basis phase
                    RMSDDegdH=RMSDDegdH+Ltemp
                    !H
                    forall(istate=1:NStates,jstate=1:NStates,istate>=jstate)
                        HdLSF_H(istate,jstate)=HdLSF_H(istate,jstate)-DegeneratePoint(ip).H(istate,jstate)
                    end forall
                    temp=dsyFrobeniusSquare(HdLSF_H,NStates)
                    RMSDDegH=RMSDDegH+temp
                    L=L+DegeneratePoint(ip).weight*(Ltemp+HdLSF_EnergyScaleSquare*temp)
                    HdLSF_dcH=sy3UnitaryTransformation(dcHd_ByKnownf(HdLSF_f),HdLSF_phi,NExpansionCoefficients,NStates)
                    do istate=1,NStates
                        HdLSF_y( indicerow:indicerow+NStates-istate)=HdLSF_EnergyScaleSquare*DegeneratePoint(ip).H(istate:NStates,istate)
                        HdLSF_MT(indicerow:indicerow+NStates-istate,:)=HdLSF_EnergyScaleSquare*transpose(HdLSF_dcH(:,istate:NStates,istate))
                        indicerow=indicerow+NStates-istate+1
                    end do
                    !▽H (▽_c phi is neglected)
                    HdLSF_dcdH=sy4UnitaryTransformation(dcdHd_ByKnownfdT(transpose(HdLSF_fd)),HdLSF_phi,NExpansionCoefficients,InternalDimension,NStates)
                    do istate=1,NStates
                        do jstate=istate,NStates
                            HdLSF_y( indicerow:indicerow+InternalDimension-1)=DegeneratePoint(ip).dH(:,jstate,istate)
                            HdLSF_MT(indicerow:indicerow+InternalDimension-1,:)=transpose(HdLSF_dcdH(:,:,jstate,istate))
                            indicerow=indicerow+InternalDimension
                        end do
                    end do
                end do
                RMSDDegH=dSqrt(RMSDDegH/(NStates*NStates)/NDegeneratePoints)
                RMSDDegdH=dSqrt(RMSDDegdH/(InternalDimension*NStates*NStates)/NDegeneratePoints)
                Ltemp=0d0
                do ip=1,NArtifactPoints!Unreliable data points, energy only
                    call AdiabaticEnergy_State_f(ArtifactPoint(ip).geom,HdLSF_energy,HdLSF_phi,HdLSF_f)
                    HdLSF_energy=HdLSF_energy-ArtifactPoint(ip).energy
                    Ltemp=Ltemp+ArtifactPoint(ip).weight*dot_product(HdLSF_energy,HdLSF_energy)
                    HdLSF_dcH=sy3UnitaryTransformation(dcHd_ByKnownf(HdLSF_f),HdLSF_phi,NExpansionCoefficients,NStates)
                    HdLSF_y(indicerow:indicerow+NStates-1)=HdLSF_EnergyScaleSquare*ArtifactPoint(ip).energy
                    forall(istate=1:NStates)
                        HdLSF_MT(indicerow+istate-1,:)=HdLSF_EnergyScaleSquare*HdLSF_dcH(:,istate,istate)
                    end forall
                    indicerow=indicerow+NStates
                end do
                L=L+HdLSF_EnergyScaleSquare*Ltemp
                !Done construction, put them into A and b
                HdLSF_M=transpose(HdLSF_MT)
                forall(ip=1:HdLSF_NData)
                    HdLSF_M(:,ip)=HdLSF_M(:,ip)*HdLSF_W(ip)
                end forall
                A=matmul(HdLSF_M,HdLSF_MT)
                forall(ip=1:NExpansionCoefficients)!Regularization
                    A(ip,ip)=A(ip,ip)+HdLSF_Regularization
                end forall
                b=matmul(HdLSF_M,HdLSF_y)
            end subroutine LSFMatrices_L_RMSD
    end subroutine pseudolinear
    
    subroutine TrustRegionInterface(c)!Fit Hd by trust region method
        !To save CPU time, weight -> Sqrt(weight), and will recover on exit
        real*8,dimension(NExpansionCoefficients),intent(inout)::c
        integer::i
        real*8::L,RMSDenergy,RMSDdH,RMSDDegH,RMSDDegdH
        !Initialize
            !weight -> Sqrt(weight)
                forall(i=1:NPoints)
                    point(i).weight=dSqrt(point(i).weight)
                end forall
                forall(i=1:NDegeneratePoints)
                    DegeneratePoint(i).weight=dSqrt(DegeneratePoint(i).weight)
                end forall
                forall(i=1:NArtifactPoints)
                    ArtifactPoint(i).weight=dSqrt(ArtifactPoint(i).weight)
                end forall
            !Allocate global local minimizer work space
                if(allocated(HdLSF_dHd)) deallocate(HdLSF_dHd)
                allocate(HdLSF_dHd(InternalDimension,NStates,NStates))
                if(allocated(HdLSF_dcphi)) deallocate(HdLSF_dcphi)
                allocate(HdLSF_dcphi(NExpansionCoefficients,NStates,NStates))
                if(allocated(HdLSF_dcHrep)) deallocate(HdLSF_dcHrep)
                allocate(HdLSF_dcHrep(NExpansionCoefficients,NStates,NStates))
                if(allocated(HdLSF_dcdHrep)) deallocate(HdLSF_dcdHrep)
                allocate(HdLSF_dcdHrep(NExpansionCoefficients,InternalDimension,NStates,NStates))
                !Trust region
                    if(allocated(HdLSF_spResidue)) deallocate(HdLSF_spResidue)
                    allocate(HdLSF_spResidue(DataPerPoint))
                    if(allocated(HdLSF_sdegpResidue)) deallocate(HdLSF_sdegpResidue)
                    allocate(HdLSF_sdegpResidue(DataPerDegeneratePoint))
                    if(allocated(HdLSF_spJacobian)) deallocate(HdLSF_spJacobian)
                    allocate(HdLSF_spJacobian(DataPerPoint,NExpansionCoefficients))
                    if(allocated(HdLSF_sdegpJacobian)) deallocate(HdLSF_sdegpJacobian)
                    allocate(HdLSF_sdegpJacobian(DataPerDegeneratePoint,NExpansionCoefficients))
        !Solve
        call showtime()
        write(*,'(1x,A40)')'Search for local minimum by trust region'
        call TrustRegion(Residue,c,HdLSF_NData+NExpansionCoefficients,NExpansionCoefficients,Jacobian=Jacobian,&
            MaxIteration=HdLSF_MaxLocalMinimizerIteration)
        !Clean up
            !weight <- Sqrt(weight)
                forall(i=1:NPoints)
                    point(i).weight=point(i).weight*point(i).weight
                end forall
                forall(i=1:NDegeneratePoints)
                    DegeneratePoint(i).weight=DegeneratePoint(i).weight*DegeneratePoint(i).weight
                end forall
                forall(i=1:NArtifactPoints)
                    ArtifactPoint(i).weight=ArtifactPoint(i).weight*ArtifactPoint(i).weight
                end forall
            !Global local minimizer work space
                deallocate(HdLSF_dHd)
                deallocate(HdLSF_dcphi)
                deallocate(HdLSF_dcHrep)
                deallocate(HdLSF_dcdHrep)
                !Trust region
                    deallocate(HdLSF_spResidue)
                    deallocate(HdLSF_sdegpResidue)
                    deallocate(HdLSF_spJacobian)
                    deallocate(HdLSF_sdegpJacobian)
        !Output
        call L_RMSD(c,L,RMSDenergy,RMSDdH,RMSDDegH,RMSDDegdH)
        write(*,'(1x,A23)')'Result of trust region:'
        write(*,*)'Lagrangian =',L
        write(*,*)'RMSD over regular data points:'
        write(*,*)'     E =',RMSDenergy/cm_1InAu,'cm-1'
        write(*,*)'    dH =',RMSDdH,'a.u.'
        if(NDegeneratePoints>0) then
            write(*,*)'RMSD over almost degenerate data points:'
            write(*,*)'     H =',RMSDDegH/cm_1InAu,'cm-1'
            write(*,*)'    dH =',RMSDDegdH,'a.u.'
        end if
        write(*,'(1x,A30)')'Save Hd expansion coefficients'
        call c2HdEC(c,HdEC,NExpansionCoefficients)
        call WriteHdExpansionCoefficients(HdEC,NStates,NOrder)
        contains
            !To save memory, off-diagonals are treated as twice weighed
            subroutine Residue(r,c,M,N)
                integer,intent(in)::M,N
                real*8,dimension(M),intent(out)::r
                real*8,dimension(N),intent(in)::c
                integer::ip,istate,jstate,indicerow,indicesp
                real*8::Ltemp
                real*8,dimension(NStates,NStates)::phi
                !Initialize
                    r(HdLSF_NData+1:M)=HdLSF_SqrtRegularization*c!Regularization
                    call c2HdEC(c,HdEC,NExpansionCoefficients)
                indicerow=1!Start from 1st row
                do ip=1,NPoints!Regular data points
                    call AdiabaticEnergy_dH(point(ip).geom,HdLSF_energy,HdLSF_dH)!Adiabatic representation
                    call dFixdHPhase(HdLSF_dH,point(ip).dH,Ltemp,InternalDimension,NStates)!Fix off-diagonals phase
                    HdLSF_spResidue(1:NStates)=HdLSF_EnergyScale*(HdLSF_energy-point(ip).energy)!Energy
                    indicesp=NStates+1
                    do istate=1,NStates!▽H
                        HdLSF_spResidue(indicesp:indicesp+InternalDimension-1)=HdLSF_dH(:,istate,istate)-point(ip).dH(:,istate,istate)!Diagonal
                        indicesp=indicesp+InternalDimension
                        do jstate=istate+1,NStates!Off-diagonal
                            HdLSF_spResidue(indicesp:indicesp+InternalDimension-1)=Sqrt2*(HdLSF_dH(:,jstate,istate)-point(ip).dH(:,jstate,istate))
                            indicesp=indicesp+InternalDimension
                        end do
                    end do
                    r(indicerow:indicerow+DataPerPoint-1)=point(ip).weight*HdLSF_spResidue
                    indicerow=indicerow+DataPerPoint
                end do
                do ip=1,NDegeneratePoints!Almost degenerate data points
                    call NondegenerateH_dH(DegeneratePoint(ip).geom,HdLSF_H,HdLSF_dH)!Nondegenerate representation
                    call dFixHPhaseBydH(HdLSF_H,HdLSF_dH,DegeneratePoint(ip).dH,Ltemp,InternalDimension,NStates)!Fix off-diagonals phase
                    indicesp=1
                    do istate=1,NStates!H
                        HdLSF_sdegpResidue(indicesp)=HdLSF_EnergyScale*(HdLSF_H(istate,istate)-DegeneratePoint(ip).H(istate,istate))!Diagonal
                        HdLSF_sdegpResidue(indicesp+1:indicesp+NStates-istate)=Sqrt2*HdLSF_EnergyScale*(HdLSF_H(istate+1:NStates,istate)-DegeneratePoint(ip).H(istate+1:NStates,istate))
                        indicesp=indicesp+NStates-istate+1
                    end do
                    do istate=1,NStates!▽H
                        HdLSF_sdegpResidue(indicesp:indicesp+InternalDimension-1)=HdLSF_dH(:,istate,istate)-DegeneratePoint(ip).dH(:,istate,istate)!Diagonal
                        indicesp=indicesp+InternalDimension
                        do jstate=istate+1,NStates!Off-diagonal
                            HdLSF_sdegpResidue(indicesp:indicesp+InternalDimension-1)=Sqrt2*(HdLSF_dH(:,jstate,istate)-DegeneratePoint(ip).dH(:,jstate,istate))
                            indicesp=indicesp+InternalDimension
                        end do
                    end do
                    r(indicerow:indicerow+DataPerDegeneratePoint-1)=DegeneratePoint(ip).weight*HdLSF_sdegpResidue
                    indicerow=indicerow+DataPerDegeneratePoint
                end do
                do ip=1,NArtifactPoints!Unreliable data points, energy only
                    r(indicerow:indicerow+NStates-1)=ArtifactPoint(ip).weight*HdLSF_EnergyScale*(AdiabaticEnergy(ArtifactPoint(ip).geom)-ArtifactPoint(ip).energy)
                    indicerow=indicerow+NStates
                end do
            end subroutine Residue
            integer function Jacobian(Jacob,c,M,N)
                integer,intent(in)::M,N
                real*8,dimension(M,N),intent(out)::Jacob
                real*8,dimension(N),intent(in)::c
                integer::ip,istate,jstate,i,indicerow,indicesp
                real*8::Ltemp
                real*8,dimension(NStates,NStates)::phi
                real*8,dimension(InternalDimension)::gradienttemp
                !Initialize
                    Jacob(HdLSF_NData+1:M,:)=0d0
                    forall(ip=1:NExpansionCoefficients)!Regularization
                        Jacob(HdLSF_NData+ip,ip)=HdLSF_SqrtRegularization
                    end forall
                    call c2HdEC(c,HdEC,NExpansionCoefficients)
                indicerow=1!Start from 1st row
                do ip=1,NPoints!Regular data points
                    call AdiabaticEnergy_dH_State_f_fd(point(ip).geom,HdLSF_energy,HdLSF_dH,HdLSF_phi,HdLSF_f,HdLSF_fd)!Adiabatic representation
                    call dAssignBasisPhaseBydH(HdLSF_phi,HdLSF_dH,point(ip).dH,Ltemp,InternalDimension,NStates)!Assign basis phase
                    !▽_c phi
                    HdLSF_dcH=sy3UnitaryTransformation(dcHd_ByKnownf(HdLSF_f),HdLSF_phi,NExpansionCoefficients,NStates)
                    HdLSF_dcphi=deigvec_ByKnowneigval_dA(HdLSF_energy,HdLSF_dcH,NExpansionCoefficients,NStates)
                    !Energy
                    forall(istate=1:NStates)
                        HdLSF_spJacobian(istate,:)=HdLSF_EnergyScale*HdLSF_dcH(:,istate,istate)
                    end forall
                    !▽H
                    HdLSF_dcdHrep=asy3matdirectmulsy3(HdLSF_dcphi,HdLSF_dH,NExpansionCoefficients,InternalDimension,NStates)
                    HdLSF_dcdHrep=sy4UnitaryTransformation(dcdHd_ByKnownfdT(transpose(HdLSF_fd)),HdLSF_phi,NExpansionCoefficients,InternalDimension,NStates)&
                        -HdLSF_dcdHrep-transpose4(HdLSF_dcdHrep,NExpansionCoefficients,InternalDimension,NStates,NStates)
                    indicesp=NStates+1
                    do istate=1,NStates
                        HdLSF_spJacobian(indicesp:indicesp+InternalDimension-1,:)=transpose(HdLSF_dcdHrep(:,:,istate,istate))!Diagonal
                        indicesp=indicesp+InternalDimension
                        do jstate=istate+1,NStates!Off-diagonal
                            HdLSF_spJacobian(indicesp:indicesp+InternalDimension-1,:)=Sqrt2*transpose(HdLSF_dcdHrep(:,:,jstate,istate))
                            indicesp=indicesp+InternalDimension
                        end do
                    end do
                    Jacob(indicerow:indicerow+DataPerPoint-1,:)=point(ip).weight*HdLSF_spJacobian
                    indicerow=indicerow+DataPerPoint
                end do
                do ip=1,NDegeneratePoints!Almost degenerate data points
                    !In this loop, HdLSF_energy stores the eigenvalues of nondegenerate operator
                    call NondegenerateH_dH_eigval_State_dHd_f_fd(DegeneratePoint(ip).geom,HdLSF_H,HdLSF_dH,HdLSF_energy,HdLSF_phi,HdLSF_dHd,HdLSF_f,HdLSF_fd)!Nondegenerate representation
                    call dFixHPhase_AssignBasisPhaseBydH(HdLSF_H,HdLSF_phi,HdLSF_dH,DegeneratePoint(ip).dH,Ltemp,InternalDimension,NStates)!Assign basis phase   
                    !▽_c phi
                    HdLSF_dcphi=deigvec_ByKnowneigval_dA(HdLSF_energy,sy3UnitaryTransformation(dcAd_ByKnown(dcdHd_ByKnownfdT(transpose(HdLSF_fd)),HdLSF_dHd),HdLSF_phi,NExpansionCoefficients,NStates),NExpansionCoefficients,NStates)
                    !H
                    HdLSF_dcHrep=asy3matmulsy(HdLSF_dcphi,HdLSF_H,NExpansionCoefficients,NStates)
                    HdLSF_dcHrep=sy3UnitaryTransformation(dcHd_ByKnownf(HdLSF_f),HdLSF_phi,NExpansionCoefficients,NStates)&
                        -HdLSF_dcHrep-transpose3(HdLSF_dcHrep,NExpansionCoefficients,NStates,NStates)
                    indicesp=1
                    do istate=1,NStates
                        HdLSF_sdegpJacobian(indicesp,:)=HdLSF_EnergyScale*HdLSF_dcHrep(:,istate,istate)!Diagonal
                        forall(jstate=istate+1:NStates)!Off-diagonal
                            HdLSF_sdegpJacobian(indicesp+jstate-istate,:)=Sqrt2*HdLSF_EnergyScale*HdLSF_dcHrep(:,jstate,istate)
                        end forall
                        indicesp=indicesp+NStates-istate+1
                    end do
                    !▽H
                    HdLSF_dcdHrep=asy3matdirectmulsy3(HdLSF_dcphi,HdLSF_dH,NExpansionCoefficients,InternalDimension,NStates)
                    HdLSF_dcdHrep=sy4UnitaryTransformation(dcdHd_ByKnownfdT(transpose(HdLSF_fd)),HdLSF_phi,NExpansionCoefficients,InternalDimension,NStates)&
                        -HdLSF_dcdHrep-transpose4(HdLSF_dcdHrep,NExpansionCoefficients,InternalDimension,NStates,NStates)
                    do istate=1,NStates
                        HdLSF_spJacobian(indicesp:indicesp+InternalDimension-1,:)=transpose(HdLSF_dcdHrep(:,:,istate,istate))!Diagonal
                        indicesp=indicesp+InternalDimension
                        do jstate=istate+1,NStates!Off-diagonal
                            HdLSF_spJacobian(indicesp:indicesp+InternalDimension-1,:)=Sqrt2*transpose(HdLSF_dcdHrep(:,:,jstate,istate))
                            indicesp=indicesp+InternalDimension
                        end do
                    end do
                    Jacob(indicerow:indicerow+DataPerDegeneratePoint-1,:)=DegeneratePoint(ip).weight*HdLSF_sdegpJacobian
                    indicerow=indicerow+DataPerDegeneratePoint
                end do
                do ip=1,NArtifactPoints!Unreliable data points, energy only
                    call AdiabaticEnergy_State_f(ArtifactPoint(ip).geom,HdLSF_energy,HdLSF_phi,HdLSF_f)
                    HdLSF_dcH=sy3UnitaryTransformation(dcHd_ByKnownf(HdLSF_f),HdLSF_phi,NExpansionCoefficients,NStates)
                    forall(istate=1:NStates)
                        Jacob(indicerow+istate-1,:)=ArtifactPoint(ip).weight*HdLSF_EnergyScale*HdLSF_dcH(:,istate,istate)
                    end forall
                    indicerow=indicerow+NStates
                end do
                Jacobian=0!return 0
            end function Jacobian
    end subroutine TrustRegionInterface
    
    subroutine LineSearchInterface(c)!Fit Hd by line search method
        real*8,dimension(NExpansionCoefficients),intent(inout)::c
        real*8::L,RMSDenergy,RMSDdH,RMSDDegH,RMSDDegdH
        !Initialize
            !Allocate global local minimizer work space
            if(allocated(HdLSF_dHd)) deallocate(HdLSF_dHd)
            allocate(HdLSF_dHd(InternalDimension,NStates,NStates))
            if(allocated(HdLSF_dcphi)) deallocate(HdLSF_dcphi)
            allocate(HdLSF_dcphi(NExpansionCoefficients,NStates,NStates))
            if(allocated(HdLSF_dcHrep)) deallocate(HdLSF_dcHrep)
            allocate(HdLSF_dcHrep(NExpansionCoefficients,NStates,NStates))
            if(allocated(HdLSF_dcdHrep)) deallocate(HdLSF_dcdHrep)
            allocate(HdLSF_dcdHrep(NExpansionCoefficients,InternalDimension,NStates,NStates))
        !Solve
        call showtime()
        select case(HdLSF_LineSearcher)
            case('LBFGS')
                write(*,'(1x,A67)')'Search for local minimum by limited memory BFGS quasi-Newton method'
                call LBFGS(Lagrangian,LagrangianGradient,c,NExpansionCoefficients,f_fd=Lagrangian_LagrangianGradient,&
                    Memory=HdLSF_LBFGSMemory,Strong=UseStrongWolfe,MaxIteration=HdLSF_MaxLocalMinimizerIteration)
            case('ConjugateGradient')
                write(*,'(1x,A53)')'Search for local minimum by conjugate gradient method'
                call ConjugateGradient(Lagrangian,LagrangianGradient,c,NExpansionCoefficients,f_fd=Lagrangian_LagrangianGradient,&
                    Method=HdLSF_ConjugateGradientSolver,Strong=UseStrongWolfe,MaxIteration=HdLSF_MaxLocalMinimizerIteration)
            case default!Throw a warning
                write(*,'(1x,A40,1x,A32)')'Program abort: unsupported line searcher',HdLSF_LineSearcher
                stop
        end select
        !Output
        call L_RMSD(c,L,RMSDenergy,RMSDdH,RMSDDegH,RMSDDegdH)
        write(*,'(1x,A15)')'Result of LBFGS:'
        write(*,*)'Lagrangian =',L
        write(*,*)'RMSD over regular data points:'
        write(*,*)'     E =',RMSDenergy/cm_1InAu,'cm-1'
        write(*,*)'    dH =',RMSDdH,'a.u.'
        if(NDegeneratePoints>0) then
            write(*,*)'RMSD over almost degenerate data points:'
            write(*,*)'     H =',RMSDDegH/cm_1InAu,'cm-1'
            write(*,*)'    dH =',RMSDDegdH,'a.u.'
        end if
        write(*,'(1x,A30)')'Save Hd expansion coefficients'
        call c2HdEC(c,HdEC,NExpansionCoefficients)
        call WriteHdExpansionCoefficients(HdEC,NStates,NOrder)
        !Clean up
            deallocate(HdLSF_dHd)
            deallocate(HdLSF_dcphi)
            deallocate(HdLSF_dcHrep)
            deallocate(HdLSF_dcdHrep)
        contains
            !To save CPU time, Lagrangian -> Lagrangian / 2
            subroutine Lagrangian(L,c,dim)!dim dimensional vector c, L harvests Lagrangian
                real*8,intent(out)::L
                integer,intent(in)::dim
                real*8,dimension(dim),intent(in)::c
                integer::ip,istate,jstate
                real*8::Ltemp
                !Initialize
                    L=HdLSF_Regularization*dot_product(c,c)!Regularization
                    call c2HdEC(c,HdEC,NExpansionCoefficients)
                do ip=1,NPoints!Regular data points
                    call AdiabaticEnergy_dH(point(ip).geom,HdLSF_energy,HdLSF_dH)!Adiabatic representation
                    call dFixdHPhase(HdLSF_dH,point(ip).dH,Ltemp,InternalDimension,NStates)!Fix off-diagonals phase
                    HdLSF_energy=HdLSF_energy-point(ip).energy!Energy (▽H has done during fixing)
                    L=L+point(ip).weight*(Ltemp+HdLSF_EnergyScaleSquare*dot_product(HdLSF_energy,HdLSF_energy))
                end do
                do ip=1,NDegeneratePoints!Almost degenerate data points
                    call NondegenerateH_dH(DegeneratePoint(ip).geom,HdLSF_H,HdLSF_dH)!Nondegenerate representation
                    call dFixHPhaseBydH(HdLSF_H,HdLSF_dH,DegeneratePoint(ip).dH,Ltemp,InternalDimension,NStates)!Fix off-diagonals phase
                    forall(istate=1:NStates,jstate=1:NStates,istate>=jstate)!H (▽H has done during fixing)
                        HdLSF_H(istate,jstate)=HdLSF_H(istate,jstate)-DegeneratePoint(ip).H(istate,jstate)
                    end forall
                    L=L+DegeneratePoint(ip).weight*(Ltemp+HdLSF_EnergyScaleSquare*dsyFrobeniusSquare(HdLSF_H,NStates))
                end do
                Ltemp=0d0
                do ip=1,NArtifactPoints!Unreliable data points, energy only
                    HdLSF_energy=AdiabaticEnergy(ArtifactPoint(ip).geom)-ArtifactPoint(ip).energy
                    Ltemp=Ltemp+ArtifactPoint(ip).weight*dot_product(HdLSF_energy,HdLSF_energy)
                end do
                L=(L+HdLSF_EnergyScaleSquare*Ltemp)/2d0
            end subroutine Lagrangian
            subroutine LagrangianGradient(Ld,c,dim)!dim dimensional vector Ld & c, Ld harvests the gradient of Lagrangian over c
                integer,intent(in)::dim
                real*8,dimension(dim),intent(out)::Ld
                real*8,dimension(dim),intent(in)::c
                integer::i,ip,istate,jstate
                real*8::Ltemp
                real*8,dimension(dim)::Ldtemp
                !Initialize
                    Ld=HdLSF_Regularization*c!Regularization
                    call c2HdEC(c,HdEC,NExpansionCoefficients)
                do ip=1,NPoints!Regular data points
                    call AdiabaticEnergy_dH_State_f_fd(point(ip).geom,HdLSF_energy,HdLSF_dH,HdLSF_phi,HdLSF_f,HdLSF_fd)!Adiabatic representation
                    call dAssignBasisPhaseBydH(HdLSF_phi,HdLSF_dH,point(ip).dH,Ltemp,InternalDimension,NStates)!Assign basis phase
                    !▽_c phi
                    HdLSF_dcH=sy3UnitaryTransformation(dcHd_ByKnownf(HdLSF_f),HdLSF_phi,NExpansionCoefficients,NStates)
                    HdLSF_dcphi=deigvec_ByKnowneigval_dA(HdLSF_energy,HdLSF_dcH,NExpansionCoefficients,NStates)
                    !Energy
                    forall(istate=1:NStates)
                        HdLSF_dcHrep(:,istate,istate)=HdLSF_dcH(:,istate,istate)
                    end forall
                    forall(istate=2:NStates,jstate=1:NStates-1,istate>jstate)
                        HdLSF_dcHrep(:,istate,jstate)=0d0
                    end forall
                    HdLSF_H=diag(HdLSF_energy-point(ip).energy,NStates)
                    !▽H
                    HdLSF_dcdHrep=asy3matdirectmulsy3(HdLSF_dcphi,HdLSF_dH,NExpansionCoefficients,InternalDimension,NStates)
                    HdLSF_dcdHrep=sy4UnitaryTransformation(dcdHd_ByKnownfdT(transpose(HdLSF_fd)),HdLSF_phi,NExpansionCoefficients,InternalDimension,NStates)&
                        -HdLSF_dcdHrep-transpose4(HdLSF_dcdHrep,NExpansionCoefficients,InternalDimension,NStates,NStates)
                    forall(istate=1:NStates,jstate=1:NStates,istate>=jstate)
                        HdLSF_dH(:,istate,jstate)=HdLSF_dH(:,istate,jstate)-point(ip).dH(:,istate,jstate)
                    end forall
                    Ld=Ld+point(ip).weight*(&
                        HdLSF_EnergyScaleSquare*Trace3(sy3matmulsy(HdLSF_dcHrep,HdLSF_H,NExpansionCoefficients,NStates),NExpansionCoefficients,NStates)&    
                        +Trace3(sy4matdotmulsy3(HdLSF_dcdHrep,HdLSF_dH,NExpansionCoefficients,InternalDimension,NStates),NExpansionCoefficients,NStates))
                end do
                do ip=1,NDegeneratePoints!Almost degenerate data points
                    !In this loop, HdLSF_energy stores the eigenvalues of nondegenerate operator
                    call NondegenerateH_dH_eigval_State_dHd_f_fd(DegeneratePoint(ip).geom,HdLSF_H,HdLSF_dH,HdLSF_energy,HdLSF_phi,HdLSF_dHd,HdLSF_f,HdLSF_fd)!Nondegenerate representation
                    call dFixHPhase_AssignBasisPhaseBydH(HdLSF_H,HdLSF_phi,HdLSF_dH,DegeneratePoint(ip).dH,Ltemp,InternalDimension,NStates)!Assign basis phase
                    !▽_c phi
                    HdLSF_dcphi=deigvec_ByKnowneigval_dA(HdLSF_energy,sy3UnitaryTransformation(dcAd_ByKnown(dcdHd_ByKnownfdT(transpose(HdLSF_fd)),HdLSF_dHd),HdLSF_phi,NExpansionCoefficients,NStates),NExpansionCoefficients,NStates)
                    !H
                    HdLSF_dcHrep=asy3matmulsy(HdLSF_dcphi,HdLSF_H,NExpansionCoefficients,NStates)
                    HdLSF_dcHrep=sy3UnitaryTransformation(dcHd_ByKnownf(HdLSF_f),HdLSF_phi,NExpansionCoefficients,NStates)&
                        -HdLSF_dcHrep-transpose3(HdLSF_dcHrep,NExpansionCoefficients,NStates,NStates)
                    forall(istate=1:NStates,jstate=1:NStates,istate>=jstate)
                        HdLSF_H(istate,jstate)=HdLSF_H(istate,jstate)-DegeneratePoint(ip).H(istate,jstate)
                    end forall
                    !▽H
                    HdLSF_dcdHrep=asy3matdirectmulsy3(HdLSF_dcphi,HdLSF_dH,NExpansionCoefficients,InternalDimension,NStates)
                    HdLSF_dcdHrep=sy4UnitaryTransformation(dcdHd_ByKnownfdT(transpose(HdLSF_fd)),HdLSF_phi,NExpansionCoefficients,InternalDimension,NStates)&
                        -HdLSF_dcdHrep-transpose4(HdLSF_dcdHrep,NExpansionCoefficients,InternalDimension,NStates,NStates)
                    forall(istate=1:NStates,jstate=1:NStates,istate>=jstate)
                        HdLSF_dH(:,istate,jstate)=HdLSF_dH(:,istate,jstate)-DegeneratePoint(ip).dH(:,istate,jstate)
                    end forall
                    Ld=Ld+DegeneratePoint(ip).weight*(&
                        HdLSF_EnergyScaleSquare*Trace3(sy3matmulsy(HdLSF_dcHrep,HdLSF_H,NExpansionCoefficients,NStates),NExpansionCoefficients,NStates)&
                        +Trace3(sy4matdotmulsy3(HdLSF_dcdHrep,HdLSF_dH,NExpansionCoefficients,InternalDimension,NStates),NExpansionCoefficients,NStates))
                end do
                Ldtemp=0d0
                do ip=1,NArtifactPoints!Unreliable data points, energy only
                    call AdiabaticEnergy_State_f(ArtifactPoint(ip).geom,HdLSF_energy,HdLSF_phi,HdLSF_f)
                    HdLSF_dcH=sy3UnitaryTransformation(dcHd_ByKnownf(HdLSF_f),HdLSF_phi,NExpansionCoefficients,NStates)
                    forall(istate=1:NStates)
                        HdLSF_dcHrep(:,istate,istate)=HdLSF_dcH(:,istate,istate)
                    end forall
                    forall(istate=2:NStates,jstate=1:NStates-1,istate>jstate)
                        HdLSF_dcHrep(:,istate,jstate)=0d0
                    end forall
                    HdLSF_H=diag(HdLSF_energy-ArtifactPoint(ip).energy,NStates)
                    Ldtemp=Ldtemp+ArtifactPoint(ip).weight*Trace3(sy3matmulsy(HdLSF_dcHrep,HdLSF_H,NExpansionCoefficients,NStates),NExpansionCoefficients,NStates)
                end do
                Ld=Ld+HdLSF_EnergyScaleSquare*Ldtemp
            end subroutine LagrangianGradient
            integer function Lagrangian_LagrangianGradient(L,Ld,c,dim)!dim dimensional vector Ld & c, L harvests Lagrangian, Ld harvests the gradient of Lagrangian over c
                integer,intent(in)::dim
                real*8,intent(out)::L
                real*8,dimension(dim),intent(out)::Ld
                real*8,dimension(dim),intent(in)::c
                integer::i,ip,istate,jstate
                real*8::Ltemp
                real*8,dimension(dim)::Ldtemp
                !Initialize
                    L =HdLSF_Regularization*dot_product(c,c)!Regularization
                    Ld=HdLSF_Regularization*c!Regularization
                    call c2HdEC(c,HdEC,NExpansionCoefficients)
                do ip=1,NPoints!Regular data points
                    call AdiabaticEnergy_dH_State_f_fd(point(ip).geom,HdLSF_energy,HdLSF_dH,HdLSF_phi,HdLSF_f,HdLSF_fd)!Adiabatic representation
                    call dAssignBasisPhaseBydH(HdLSF_phi,HdLSF_dH,point(ip).dH,Ltemp,InternalDimension,NStates)!Assign basis phase
                    !▽_c phi
                    HdLSF_dcH=sy3UnitaryTransformation(dcHd_ByKnownf(HdLSF_f),HdLSF_phi,NExpansionCoefficients,NStates)
                    HdLSF_dcphi=deigvec_ByKnowneigval_dA(HdLSF_energy,HdLSF_dcH,NExpansionCoefficients,NStates)
                    !Energy
                    forall(istate=1:NStates)
                        HdLSF_dcHrep(:,istate,istate)=HdLSF_dcH(:,istate,istate)
                    end forall
                    forall(istate=2:NStates,jstate=1:NStates-1,istate>jstate)
                        HdLSF_dcHrep(:,istate,jstate)=0d0
                    end forall
                    HdLSF_H=diag(HdLSF_energy-point(ip).energy,NStates)
                    !▽H
                    HdLSF_dcdHrep=asy3matdirectmulsy3(HdLSF_dcphi,HdLSF_dH,NExpansionCoefficients,InternalDimension,NStates)
                    HdLSF_dcdHrep=sy4UnitaryTransformation(dcdHd_ByKnownfdT(transpose(HdLSF_fd)),HdLSF_phi,NExpansionCoefficients,InternalDimension,NStates)&
                        -HdLSF_dcdHrep-transpose4(HdLSF_dcdHrep,NExpansionCoefficients,InternalDimension,NStates,NStates)
                    forall(istate=1:NStates,jstate=1:NStates,istate>=jstate)
                        HdLSF_dH(:,istate,jstate)=HdLSF_dH(:,istate,jstate)-point(ip).dH(:,istate,jstate)
                    end forall
                    L =L +point(ip).weight*(Ltemp+HdLSF_EnergyScaleSquare*dsyFrobeniusSquare(HdLSF_H,NStates))
                    Ld=Ld+point(ip).weight*(&
                        HdLSF_EnergyScaleSquare*Trace3(sy3matmulsy(HdLSF_dcHrep,HdLSF_H,NExpansionCoefficients,NStates),NExpansionCoefficients,NStates)&    
                        +Trace3(sy4matdotmulsy3(HdLSF_dcdHrep,HdLSF_dH,NExpansionCoefficients,InternalDimension,NStates),NExpansionCoefficients,NStates))
                end do
                do ip=1,NDegeneratePoints!Almost degenerate data points
                    !In this loop, HdLSF_energy stores the eigenvalues of nondegenerate operator
                    call NondegenerateH_dH_eigval_State_dHd_f_fd(DegeneratePoint(ip).geom,HdLSF_H,HdLSF_dH,HdLSF_energy,HdLSF_phi,HdLSF_dHd,HdLSF_f,HdLSF_fd)!Nondegenerate representation
                    call dFixHPhase_AssignBasisPhaseBydH(HdLSF_H,HdLSF_phi,HdLSF_dH,DegeneratePoint(ip).dH,Ltemp,InternalDimension,NStates)!Assign basis phase
                    !▽_c phi
                    HdLSF_dcphi=deigvec_ByKnowneigval_dA(HdLSF_energy,sy3UnitaryTransformation(dcAd_ByKnown(dcdHd_ByKnownfdT(transpose(HdLSF_fd)),HdLSF_dHd),HdLSF_phi,NExpansionCoefficients,NStates),NExpansionCoefficients,NStates)
                    !H
                    HdLSF_dcHrep=asy3matmulsy(HdLSF_dcphi,HdLSF_H,NExpansionCoefficients,NStates)
                    HdLSF_dcHrep=sy3UnitaryTransformation(dcHd_ByKnownf(HdLSF_f),HdLSF_phi,NExpansionCoefficients,NStates)&
                        -HdLSF_dcHrep-transpose3(HdLSF_dcHrep,NExpansionCoefficients,NStates,NStates)
                    forall(istate=1:NStates,jstate=1:NStates,istate>=jstate)
                        HdLSF_H(istate,jstate)=HdLSF_H(istate,jstate)-DegeneratePoint(ip).H(istate,jstate)
                    end forall
                    !▽H
                    HdLSF_dcdHrep=asy3matdirectmulsy3(HdLSF_dcphi,HdLSF_dH,NExpansionCoefficients,InternalDimension,NStates)
                    HdLSF_dcdHrep=sy4UnitaryTransformation(dcdHd_ByKnownfdT(transpose(HdLSF_fd)),HdLSF_phi,NExpansionCoefficients,InternalDimension,NStates)&
                        -HdLSF_dcdHrep-transpose4(HdLSF_dcdHrep,NExpansionCoefficients,InternalDimension,NStates,NStates)
                    forall(istate=1:NStates,jstate=1:NStates,istate>=jstate)
                        HdLSF_dH(:,istate,jstate)=HdLSF_dH(:,istate,jstate)-DegeneratePoint(ip).dH(:,istate,jstate)
                    end forall
                    L =L +DegeneratePoint(ip).weight*(Ltemp+HdLSF_EnergyScaleSquare*dsyFrobeniusSquare(HdLSF_H,NStates))
                    Ld=Ld+DegeneratePoint(ip).weight*(&
                        HdLSF_EnergyScaleSquare*Trace3(sy3matmulsy(HdLSF_dcHrep,HdLSF_H,NExpansionCoefficients,NStates),NExpansionCoefficients,NStates)&
                        +Trace3(sy4matdotmulsy3(HdLSF_dcdHrep,HdLSF_dH,NExpansionCoefficients,InternalDimension,NStates),NExpansionCoefficients,NStates))
                end do
                Ltemp=0d0
                Ldtemp=0d0
                do ip=1,NArtifactPoints!Unreliable data points, energy only
                    call AdiabaticEnergy_State_f(ArtifactPoint(ip).geom,HdLSF_energy,HdLSF_phi,HdLSF_f)
                    HdLSF_dcH=sy3UnitaryTransformation(dcHd_ByKnownf(HdLSF_f),HdLSF_phi,NExpansionCoefficients,NStates)
                    forall(istate=1:NStates)
                        HdLSF_dcHrep(:,istate,istate)=HdLSF_dcH(:,istate,istate)
                    end forall
                    forall(istate=2:NStates,jstate=1:NStates-1,istate>jstate)
                        HdLSF_dcHrep(:,istate,jstate)=0d0
                    end forall
                    HdLSF_H=diag(HdLSF_energy-ArtifactPoint(ip).energy,NStates)
                    Ltemp = Ltemp+ArtifactPoint(ip).weight*dsyFrobeniusSquare(HdLSF_H,NStates)
                    Ldtemp=Ldtemp+ArtifactPoint(ip).weight*Trace3(sy3matmulsy(HdLSF_dcHrep,HdLSF_H,NExpansionCoefficients,NStates),NExpansionCoefficients,NStates)
                end do
                L =(L+HdLSF_EnergyScaleSquare* Ltemp)/2d0
                Ld=Ld+HdLSF_EnergyScaleSquare*Ldtemp
                Lagrangian_LagrangianGradient=0!return 0
            end function Lagrangian_LagrangianGradient
    end subroutine LineSearchInterface
!----------------- End -----------------

!---------- Auxiliary routine ----------
    !Convert c to DiabaticHamiltonian module form of Hd expansion coefficient
    subroutine c2HdEC(c,HdEC,N)
        integer,intent(in)::N
        real*8,dimension(N),intent(in)::c
        type(HdExpansionCoefficient),dimension(NStates,NStates),intent(inout)::HdEC
        integer::i,j,istate,jstate,iorder,indice
        indice=1
        do istate=1,NStates
            do jstate=istate,NStates
                do iorder=0,NOrder
                    j=size(HdEC(jstate,istate).Order(iorder).Array)
                    HdEC(jstate,istate).Order(iorder).Array=c(indice:indice+j-1)
                    indice=indice+j
                end do
            end do
        end do
    end subroutine c2HdEC
    !Inverse conversion
    subroutine HdEC2c(HdEC,c,N)
        integer,intent(in)::N
        type(HdExpansionCoefficient),dimension(NStates,NStates),intent(in)::HdEC
        real*8,dimension(N),intent(out)::c
        integer::i,j,istate,jstate,iorder,indice
        indice=1
        do istate=1,NStates
            do jstate=istate,NStates
                do iorder=0,NOrder
                    j=size(HdEC(jstate,istate).Order(iorder).Array)
                    c(indice:indice+j-1)=HdEC(jstate,istate).Order(iorder).Array
                    indice=indice+j
                end do
            end do
        end do
    end subroutine HdEC2c

    !Lagrangian and root mean square deviations are the evaluation standards for the fit
    !Input:  current c
    !Output: L harvests Lagrangian,
    !        RMSDenergy/dH harvests root mean square deviation of adiabatic energy/dH over point,
    !        RMSDDegH/dH harvests root mean square deviation of nondegenerate H/dH over DegeneratePoint
    subroutine L_RMSD(c,L,RMSDenergy,RMSDdH,RMSDDegH,RMSDDegdH)
        real*8,dimension(NExpansionCoefficients),intent(in)::c
        real*8,intent(out)::L,RMSDenergy,RMSDdH,RMSDDegH,RMSDDegdH
        integer::ip,istate,jstate
        real*8::Ltemp,temp
        !Initialize
            L=HdLSF_Regularization*dot_product(c,c)!Regularization
            call c2HdEC(c,HdEC,NExpansionCoefficients)
        RMSDenergy=0d0
        RMSDdH=0d0
        do ip=1,NPoints!Regular data points, compute RMSD
            call AdiabaticEnergy_dH(point(ip).geom,HdLSF_energy,HdLSF_dH)!Adiabatic representation
            call dFixdHPhase(HdLSF_dH,point(ip).dH,Ltemp,InternalDimension,NStates)!Fix off-diagonals phase
            RMSDdH=RMSDdH+Ltemp
            HdLSF_energy=HdLSF_energy-point(ip).energy!Energy (dH has done during fixing)
            temp=dot_product(HdLSF_energy,HdLSF_energy)
            RMSDenergy=RMSDenergy+temp
            L=L+point(ip).weight*(Ltemp+HdLSF_EnergyScaleSquare*temp)
        end do
        RMSDenergy=dSqrt(RMSDenergy/NStates/NPoints)
        RMSDdH=dSqrt(RMSDdH/(InternalDimension*NStates*NStates)/NPoints)
        RMSDDegH=0d0
        RMSDDegdH=0d0
        do ip=1,NDegeneratePoints!Almost degenerate data points, compute RMSDDeg
            call NondegenerateH_dH(DegeneratePoint(ip).geom,HdLSF_H,HdLSF_dH)!Nondegenerate representation
            call dFixHPhaseBydH(HdLSF_H,HdLSF_dH,DegeneratePoint(ip).dH,Ltemp,InternalDimension,NStates)!Fix off-diagonals phase
            RMSDDegdH=RMSDDegdH+Ltemp
            forall(istate=1:NStates,jstate=1:NStates,istate>=jstate)!H (dH has done during fixing)
                HdLSF_H(istate,jstate)=HdLSF_H(istate,jstate)-DegeneratePoint(ip).H(istate,jstate)
            end forall
            temp=dsyFrobeniusSquare(HdLSF_H,NStates)
            RMSDDegH=RMSDDegH+temp
            L=L+DegeneratePoint(ip).weight*(Ltemp+HdLSF_EnergyScaleSquare*temp)
        end do
        RMSDDegH=dSqrt(RMSDDegH/(NStates*NStates)/NDegeneratePoints)
        RMSDDegdH=dSqrt(RMSDDegdH/(InternalDimension*NStates*NStates)/NDegeneratePoints)
        Ltemp=0d0
        do ip=1,NArtifactPoints!Unreliable data points, energy only
            HdLSF_energy=AdiabaticEnergy(ArtifactPoint(ip).geom)-ArtifactPoint(ip).energy
            Ltemp=Ltemp+ArtifactPoint(ip).weight*dot_product(HdLSF_energy,HdLSF_energy)
        end do
        L=L+HdLSF_EnergyScaleSquare*Ltemp
    end subroutine L_RMSD

    !The value of ▽_cHd in diabatic representation at some coordinate q, where c is the expansion coefficient vector
    !f stores expansion basis function values at this q
    function dcHd_ByKnownf(f)
        real*8,dimension(NExpansionCoefficients,NStates,NStates)::dcHd_ByKnownf
        real*8,dimension(NExpansionBasis),intent(in)::f
        integer::i,j,indice
        indice=1
        do j=1,NStates
            do i=j,NStates
                dcHd_ByKnownf(1:indice-1,i,j)=0d0
                dcHd_ByKnownf(indice:indice+NExpansionBasis-1,i,j)=f
                indice=indice+NExpansionBasis
                dcHd_ByKnownf(indice:NExpansionCoefficients,i,j)=0d0
            end do
        end do
    end function dcHd_ByKnownf

    !The value of ▽_c▽Hd in diabatic representation at some coordinate q, where c is the expansion coefficient vector
    !fdT(i,:) stores the gradient of i-th expansion basis function
    function dcdHd_ByKnownfdT(fdT)
        real*8,dimension(NExpansionCoefficients,InternalDimension,NStates,NStates)::dcdHd_ByKnownfdT
        real*8,dimension(NExpansionBasis,InternalDimension),intent(in)::fdT
        integer::i,j,indice
        indice=1
        do j=1,NStates
            do i=j,NStates
                dcdHd_ByKnownfdT(1:indice-1,:,i,j)=0d0
                dcdHd_ByKnownfdT(indice:indice+NExpansionBasis-1,:,i,j)=fdT
                indice=indice+NExpansionBasis
                dcdHd_ByKnownfdT(indice:NExpansionCoefficients,:,i,j)=0d0
            end do
        end do
    end function dcdHd_ByKnownfdT

    !The value of ▽_cA in diabatic representation at some coordinate q, where c is the expansion coefficient vector
    !For A adopted here, i.e. (▽H)^2, we can use known ▽_c▽Hd & ▽Hd to calculate
    function dcAd_ByKnown(dcdHd,dHd)
        real*8,dimension(NExpansionCoefficients,NStates,NStates)::dcAd_ByKnown
        real*8,dimension(NExpansionCoefficients,InternalDimension,NStates,NStates),intent(in)::dcdHd
        real*8,dimension(InternalDimension,NStates,NStates),intent(in)::dHd
        dcAd_ByKnown=sy4matdotmulsy3(dcdHd,dHd,NExpansionCoefficients,InternalDimension,NStates)
        dcAd_ByKnown=dcAd_ByKnown+transpose3(dcAd_ByKnown,NExpansionCoefficients,NStates,NStates)
    end function dcAd_ByKnown
!----------------- End -----------------

end module HdLeastSquareFit