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
!where rho is the unit converter from energy to energy gradient (LSF_EnergyScale),
!    subscript ad means adiabatic representation, nd means nondegenerate representation, F means Frobenius norm
!    superscript d means diabatz, ab means ab initio
!    definition of point, DegeneratePoint, ArtifactPoint see module Basic
!We may also add a regularization to Lagrangian, namely tikhonov regularization: tau * || c ||_2^2
!where tau is the parameterized KKT multiplier (LSF_Regularization), c is undetermined parameter vector in Hd
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
        real*8::LSF_Regularization=0d0!Instead of solving the KKT multiplier, let it be a parameter
        !Available solvers: pseudolinear, TrustRegion, LBFGS, ConjugateGradient
        !Choose the nonlinear optimization solver: a single solver or a 2-step solver
        !    pseudolinear_X combining pseudolinear and another solver X
        !If you see insufficient memory, use only LBFGS or ConjugateGradient
        !Warning: pseudolinear is a quick hopper by vanishing the linear part of the gradient
        !    However, it is not necessarily able to solve the fitting alone, 
        !    unless the minimum coincidentally also has the nonlinear part of the gradient = 0
        character*32::LSF_Solver='pseudolinear_TrustRegion'
        !Max ineration control. Hopper = pseudolinear. LocalMinimizer = TrustRegion, LBFGS, ConjugateGradient
        integer::LSF_MaxHopperIteration=100,LSF_MaxLocalMinimizerIteration=1000,LSF_Max2StepIteration=10
    !pseudolinear:
        real*8::LSF_pseudolinearTol=1d-30!Convergence standard: || c_new - c_old ||^2 < absolute tolerance
        integer::LSF_pseudolinearFollowFreq=1!Every how many steps print fitting progress
    !LBFGS:
        !Choose a specific solver: LBFGS, LBFGS_Strong
        character*32::LSF_LBFGSSolver='LBFGS_Strong'
        integer::LSF_LBFGSMemory=10!Memory usage control. [ 3, 30 ] is recommended
    !ConjugateGradient:
        !Choose a specific solver: DY, DY_Strong, PR
        character*32::LSF_CGSolver='PR'

!HdLeastSquareFit module only variable
    integer::LSF_NData!Number of fitting data
    real*8::LSF_EnergyScale!Scale energy residue in Lagrangian for its unit difference from gradient
    real*8::LSF_SqrtRegularization,LSF_EnergyScaleSquare
    !Work space
        real*8,allocatable,dimension(:)::LSF_energy
        !dc = ▽_c. H, dH, dcH, dcdH are in representation of the data point 
        !For dcH & dcdH, ▽_c only operates on H & ▽H
        !phi is the basis matrix of the data point representation in diabatic representation
        real*8,allocatable,dimension(:,:)::LSF_H,LSF_phi
        real*8,allocatable,dimension(:,:,:)::LSF_dH,LSF_dcH
        real*8,allocatable,dimension(:,:,:,:)::LSF_dcdH
        real*8,allocatable,dimension(:)::LSF_f!Stores expansion basis function values
        real*8,allocatable,dimension(:,:)::LSF_fd!Stores expansion basis function gradient values
        !Pseudolinear: W is diagonal so only store its diagonal vector, MT = M^T, NT = N^T
            real*8,allocatable,dimension(:)::LSF_y,LSF_W
            real*8,allocatable,dimension(:,:)::LSF_M,LSF_MT
        !Local minimizer
            !dHd is in diabatic representaton
            !For description of dcphi, see Nonadiabatic.deigvec_ByKnowneigval_dA (with ▽ replaced by ▽_c)
            !For dcHrep & dcdHrep, ▽_c also operates on the representation basis
            real*8,allocatable,dimension(:,:,:)::LSF_dHd,LSF_dcphi,LSF_dcHrep
            real*8,allocatable,dimension(:,:,:,:)::LSF_dcdHrep
            !Trust region
                real*8,allocatable,dimension(:)::LSF_spResidue,LSF_sdegpResidue!Residue of a single (degenerate) data point
                real*8,allocatable,dimension(:,:)::LSF_spJacobian,LSF_sdegpJacobian!Jacobian of a single (degenerate) data point

contains
!The initializer for HdLeastSquareFit module
subroutine InitializeHdLeastSquareFit()
    integer::ip,istate,jstate
    real*8::MaxEnergy,MaxGrad,temp
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
    LSF_EnergyScale=MaxGrad/MaxEnergy
    LSF_EnergyScaleSquare=LSF_EnergyScale*LSF_EnergyScale
    LSF_SqrtRegularization=Sqrt(LSF_Regularization)
end subroutine InitializeHdLeastSquareFit

!Fit Hd with the designated solver
subroutine FitHd()
    integer::istate,jstate
    real*8,allocatable,dimension(:)::c
    !Initialize
        LSF_NData=DataPerPoint*NPoints+DataPerDegeneratePoint*NDegeneratePoints+NStates*NArtifactPoints
        !Initial value of c
        allocate(c(NExpansionCoefficients))
        call HdEC2c(HdEC,c,NExpansionCoefficients)
        !Allocate global work space
            if(allocated(LSF_energy)) deallocate(LSF_energy)
            allocate(LSF_energy(NStates))
            if(allocated(LSF_H)) deallocate(LSF_H)
            allocate(LSF_H(NStates,NStates))
            if(allocated(LSF_phi)) deallocate(LSF_phi)
            allocate(LSF_phi(NStates,NStates))
            if(allocated(LSF_dH)) deallocate(LSF_dH)
            allocate(LSF_dH(InternalDimension,NStates,NStates))
            if(allocated(LSF_dcH)) deallocate(LSF_dcH)
            allocate(LSF_dcH(NExpansionCoefficients,NStates,NStates))
            if(allocated(LSF_dcdH)) deallocate(LSF_dcdH)
            allocate(LSF_dcdH(NExpansionCoefficients,InternalDimension,NStates,NStates))
            if(allocated(LSF_f)) deallocate(LSF_f)
            allocate(LSF_f(NExpansionBasis))
            if(allocated(LSF_fd)) deallocate(LSF_fd)
            allocate(LSF_fd(InternalDimension,NExpansionBasis))
    !Fit Hd
    select case(LSF_Solver)
        !Single solvers
        case('pseudolinear')
            call pseudolinear(c)
        case('TrustRegion')
            call TrustRegion(c)
        case('LBFGS')
            call LimitedMemoryBFGS(c)
        case('ConjugateGradient')
            call ConjugateGradient(c)
        !2-step solvers
        case('pseudolinear_TrustRegion')
            do istate=1,LSF_Max2StepIteration
                call pseudolinear(c)
                call TrustRegion(c)
            end do
        case('pseudolinear_LBFGS')
            do istate=1,LSF_Max2StepIteration
                call pseudolinear(c)
                call LimitedMemoryBFGS(c)
            end do
        case('pseudolinear_ConjugateGradient')
            do istate=1,LSF_Max2StepIteration
                call pseudolinear(c)
                call ConjugateGradient(c)
            end do
        case default!Throw a warning
            write(*,'(1x,A50,1x,A32)')'Program abort: unsupported least square fit solver',LSF_Solver
            stop
    end select
    !Clean up
        deallocate(c)
        !Global work space
            deallocate(LSF_energy)
            deallocate(LSF_H)
            deallocate(LSF_phi)
            deallocate(LSF_dH)
            deallocate(LSF_dcH)
            deallocate(LSF_dcdH)
            deallocate(LSF_f)
            deallocate(LSF_fd)
end subroutine FitHd

!------------ Pseudolinear hopper and dependency -------------
    !To save memory, in this section off-diagonals are treated as twice weighed
    !
    !Nomenclature:
    !The general form of weighted linear least square fit with pseudoregularization is:
    ! A c = b, where A = M . W . M^T + tau, b = M . W . y,
    ! c & tau have been explained at header, y is the data vector,
    ! M^T . c is the fitting prediction of y, W is the weight

    !Fit Hd by pseudolinear method
    subroutine pseudolinear(cmin)
        real*8,dimension(NExpansionCoefficients),intent(inout)::cmin
        integer::indice,ip,i
        real*8::cchange,L,Lmin,RMSDenergy,RMSDdH,RMSDDegH,RMSDDegdH
        real*8,allocatable,dimension(:)::b,c
        real*8,allocatable,dimension(:,:)::A
        !Initialize
            !Allocate global pseudolinear work space
                if(allocated(LSF_y)) deallocate(LSF_y)
                allocate(LSF_y(LSF_NData))
                if(allocated(LSF_W)) deallocate(LSF_W)
                allocate(LSF_W(LSF_NData))
                !W will not change through out the solving procedure, so fill in its value now
                    indice=1
                    do ip=1,NPoints
                        LSF_W(indice:indice+NStates-1)=point(ip).weight!Energy
                        indice=indice+NStates
                        do i=1,NStates!▽H
                            LSF_W(indice:indice+InternalDimension-1)=point(ip).weight!Diagonal
                            indice=indice+InternalDimension
                            LSF_W(indice:indice+InternalDimension*(NStates-i)-1)=2d0*point(ip).weight!Off-diagonal
                            indice=indice+InternalDimension*(NStates-i)
                        end do
                    end do
                    do ip=1,NDegeneratePoints
                        do i=1,NStates!H
                            LSF_W(indice)=DegeneratePoint(ip).weight!Diagonal
                            LSF_W(indice+1:indice+NStates-i)=2d0*DegeneratePoint(ip).weight!Off-diagonal
                            indice=indice+NStates-i+1
                        end do
                        do i=1,NStates!▽H
                            LSF_W(indice:indice+InternalDimension-1)=DegeneratePoint(ip).weight!Diagonal
                            indice=indice+InternalDimension
                            LSF_W(indice:indice+InternalDimension*(NStates-i)-1)=2d0*DegeneratePoint(ip).weight!Off-diagonal
                            indice=indice+InternalDimension*(NStates-i)
                        end do
                    end do
                    do ip=1,NArtifactPoints
                        LSF_W(indice:indice+NStates-1)=ArtifactPoint(ip).weight!Energy only
                        indice=indice+NStates
                    end do
                if(allocated(LSF_M)) deallocate(LSF_M)
                allocate(LSF_M(NExpansionCoefficients,LSF_NData))
                if(allocated(LSF_MT)) deallocate(LSF_MT)
                allocate(LSF_MT(LSF_NData,NExpansionCoefficients))
            !Initialize linear least square fit and Lagrangian minimum
                allocate(c(NExpansionCoefficients))
                c=cmin
                allocate(b(NExpansionCoefficients))
                b=c
                allocate(A(NExpansionCoefficients,NExpansionCoefficients))
CALL L_RMSD(CMIN,L,RMSDENERGY,RMSDDH,RMSDDEGH,RMSDDEGDH)
WRITE(*,*)'L_RMSD',L
CALL LAGRANGIAN(L,CMIN,NExpansionCoefficients)
WRITE(*,*)'LAGRANGIAN',L*2D0
CALL LSFMatrices_L_RMSD(A,B,L,RMSDENERGY,RMSDDH,RMSDDEGH,RMSDDEGDH)
WRITE(*,*)'LSFMatrices_L_RMSD',L
B=C
                call LSFMatrices_L(A,b,Lmin)
WRITE(*,*)'LSFMatrices_L',LMIN
        !Solve
        call showtime()
        write(*,'(1x,A43)')'Explore phase space by pseudolinear hopping'
        do i=1,LSF_MaxHopperIteration
            call My_dposv(A,b,NExpansionCoefficients)
            cchange=dot_product(b-c,b-c)
            if(cchange<LSF_pseudolinearTol) then
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
            else if(mod(i,LSF_pseudolinearFollowFreq)==0) then
                c=b
                call LSFMatrices_L_RMSD(A,b,L,RMSDenergy,RMSDdH,RMSDDegH,RMSDDegdH)
                if(L<Lmin) then
                    Lmin=L
                    cmin=c
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
            end if
        end do
        !Clean up
            !Local work space
                deallocate(b)
                deallocate(c)
                deallocate(A)
            !Global pseudolinear work space
                deallocate(LSF_y)
                deallocate(LSF_W)
                deallocate(LSF_M)
                deallocate(LSF_MT)
        !Output
        if(i>LSF_MaxHopperIteration) then!Throw a warning
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
    end subroutine pseudolinear

    !Input:  b = current c
    !Output: A harvests A, b harvests b, L harvests Lagrangian
    subroutine LSFMatrices_L(A,b,L)
        real*8,dimension(NExpansionCoefficients,NExpansionCoefficients),intent(out)::A
        real*8,dimension(NExpansionCoefficients),intent(inout)::b
        real*8,intent(out)::L
        integer::ip,istate,jstate,indicerow
        real*8::Ltemp
        !Initialize
            L=LSF_Regularization*dot_product(b,b)!Regularization
            call c2HdEC(b,HdEC,NExpansionCoefficients)
        !Construct M^T and y, add least square fit penalty to Lagrangian
        indicerow=1!Start from 1st row
        do ip=1,NPoints!Regular data points
            call AdiabaticEnergy_dH_State_f_fd(point(ip).geom,LSF_energy,LSF_dH,LSF_phi,LSF_f,LSF_fd)!Adiabatic representation
            call dAssignBasisPhaseBydH(LSF_phi,LSF_dH,point(ip).dH,Ltemp,InternalDimension,NStates)!Assign basis phase
            !Energy
            LSF_energy=LSF_energy-point(ip).energy
            L=L+point(ip).weight*(Ltemp+LSF_EnergyScaleSquare*dot_product(LSF_energy,LSF_energy))
            LSF_dcH=sy3UnitaryTransformation(dcHd_ByKnownf(LSF_f),LSF_phi,NExpansionCoefficients,NStates)
            LSF_y(indicerow:indicerow+NStates-1)=LSF_EnergyScaleSquare*point(ip).energy
            forall(istate=1:NStates)
                LSF_MT(indicerow+istate-1,:)=LSF_EnergyScaleSquare*LSF_dcH(:,istate,istate)
            end forall
            indicerow=indicerow+NStates
            !▽H (▽_c phi is neglected)
            LSF_dcdH=sy4UnitaryTransformation(dcdHd_ByKnownfdT(transpose(LSF_fd)),LSF_phi,NExpansionCoefficients,InternalDimension,NStates)
            do istate=1,NStates
                do jstate=istate,NStates
                    LSF_y( indicerow:indicerow+InternalDimension-1)=point(ip).dH(:,jstate,istate)
                    LSF_MT(indicerow:indicerow+InternalDimension-1,:)=transpose(LSF_dcdH(:,:,jstate,istate))
                    indicerow=indicerow+InternalDimension
                end do
            end do
        end do
        do ip=1,NDegeneratePoints!Almost degenerate data points
            call NondegenerateH_dH_State_f_fd(DegeneratePoint(ip).geom,LSF_H,LSF_dH,LSF_phi,LSF_f,LSF_fd)!Nondegenerate representation
            call dFixHPhase_AssignBasisPhaseBydH(LSF_H,LSF_phi,LSF_dH,DegeneratePoint(ip).dH,Ltemp,InternalDimension,NStates)!Assign basis phase
            !H
            forall(istate=1:NStates,jstate=1:NStates,istate>=jstate)
                LSF_H(istate,jstate)=LSF_H(istate,jstate)-DegeneratePoint(ip).H(istate,jstate)
            end forall
            L=L+DegeneratePoint(ip).weight*(Ltemp+LSF_EnergyScaleSquare*dsyFrobeniusSquare(LSF_H,NStates))
            LSF_dcH=sy3UnitaryTransformation(dcHd_ByKnownf(LSF_f),LSF_phi,NExpansionCoefficients,NStates)
            do istate=1,NStates
                LSF_y( indicerow:indicerow+NStates-istate)=LSF_EnergyScaleSquare*DegeneratePoint(ip).H(istate:NStates,istate)
                LSF_MT(indicerow:indicerow+NStates-istate,:)=LSF_EnergyScaleSquare*transpose(LSF_dcH(:,istate:NStates,istate))
                indicerow=indicerow+NStates-istate+1
            end do
            !▽H (▽_c phi is neglected)
            LSF_dcdH=sy4UnitaryTransformation(dcdHd_ByKnownfdT(transpose(LSF_fd)),LSF_phi,NExpansionCoefficients,InternalDimension,NStates)
            do istate=1,NStates
                do jstate=istate,NStates
                    LSF_y( indicerow:indicerow+InternalDimension-1)=DegeneratePoint(ip).dH(:,jstate,istate)
                    LSF_MT(indicerow:indicerow+InternalDimension-1,:)=transpose(LSF_dcdH(:,:,jstate,istate))
                    indicerow=indicerow+InternalDimension
                end do
            end do
        end do
        Ltemp=0d0
        do ip=1,NArtifactPoints!Unreliable data points, energy only
            call AdiabaticEnergy_State_f(ArtifactPoint(ip).geom,LSF_energy,LSF_phi,LSF_f)
            LSF_energy=LSF_energy-ArtifactPoint(ip).energy
            Ltemp=Ltemp+ArtifactPoint(ip).weight*dot_product(LSF_energy,LSF_energy)
            LSF_dcH=sy3UnitaryTransformation(dcHd_ByKnownf(LSF_f),LSF_phi,NExpansionCoefficients,NStates)
            LSF_y(indicerow:indicerow+NStates-1)=LSF_EnergyScaleSquare*ArtifactPoint(ip).energy
            forall(istate=1:NStates)
                LSF_MT(indicerow+istate-1,:)=LSF_EnergyScaleSquare*LSF_dcH(:,istate,istate)
            end forall
            indicerow=indicerow+NStates
        end do
        L=L+LSF_EnergyScaleSquare*Ltemp
        !Done construction, put them into A and b
        LSF_M=transpose(LSF_MT)
        forall(ip=1:LSF_NData)
            LSF_M(:,ip)=LSF_M(:,ip)*LSF_W(ip)
        end forall
        A=matmul(LSF_M,LSF_MT)
        forall(ip=1:NExpansionCoefficients)!Regularization
            A(ip,ip)=A(ip,ip)+LSF_Regularization
        end forall
        b=matmul(LSF_M,LSF_y)
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
            L=LSF_Regularization*dot_product(b,b)!Regularization
            call c2HdEC(b,HdEC,NExpansionCoefficients)
        !Construct M^T and y, add least square fit penalty to Lagrangian
        indicerow=1!Start from 1st row
        RMSDenergy=0d0
        RMSDdH=0d0
        do ip=1,NPoints!Regular data points, compute RMSD
            call AdiabaticEnergy_dH_State_f_fd(point(ip).geom,LSF_energy,LSF_dH,LSF_phi,LSF_f,LSF_fd)!Adiabatic representation
            call dAssignBasisPhaseBydH(LSF_phi,LSF_dH,point(ip).dH,Ltemp,InternalDimension,NStates)!Assign basis phase
            RMSDdH=RMSDdH+Ltemp
            !Energy
            LSF_energy=LSF_energy-point(ip).energy
            temp=dot_product(LSF_energy,LSF_energy)
            RMSDenergy=RMSDenergy+temp
            L=L+point(ip).weight*(Ltemp+LSF_EnergyScaleSquare*temp)
            LSF_dcH=sy3UnitaryTransformation(dcHd_ByKnownf(LSF_f),LSF_phi,NExpansionCoefficients,NStates)
            LSF_y(indicerow:indicerow+NStates-1)=LSF_EnergyScaleSquare*point(ip).energy
            forall(istate=1:NStates)
                LSF_MT(indicerow+istate-1,:)=LSF_EnergyScaleSquare*LSF_dcH(:,istate,istate)
            end forall
            indicerow=indicerow+NStates
            !▽H (▽_c phi is neglected)
            LSF_dcdH=sy4UnitaryTransformation(dcdHd_ByKnownfdT(transpose(LSF_fd)),LSF_phi,NExpansionCoefficients,InternalDimension,NStates)
            do istate=1,NStates
                do jstate=istate,NStates
                    LSF_y( indicerow:indicerow+InternalDimension-1)=point(ip).dH(:,jstate,istate)
                    LSF_MT(indicerow:indicerow+InternalDimension-1,:)=transpose(LSF_dcdH(:,:,jstate,istate))
                    indicerow=indicerow+InternalDimension
                end do
            end do
        end do
        RMSDenergy=Sqrt(RMSDenergy/NStates/NPoints)
        RMSDdH=Sqrt(RMSDdH/(InternalDimension*NStates*NStates)/NPoints)
        RMSDDegH=0d0
        RMSDDegdH=0d0
        do ip=1,NDegeneratePoints!Almost degenerate data points, compute RMSDDeg
            call NondegenerateH_dH_State_f_fd(DegeneratePoint(ip).geom,LSF_H,LSF_dH,LSF_phi,LSF_f,LSF_fd)!Nondegenerate representation
            call dFixHPhase_AssignBasisPhaseBydH(LSF_H,LSF_phi,LSF_dH,DegeneratePoint(ip).dH,Ltemp,InternalDimension,NStates)!Assign basis phase
            RMSDDegdH=RMSDDegdH+Ltemp
            !H
            forall(istate=1:NStates,jstate=1:NStates,istate>=jstate)
                LSF_H(istate,jstate)=LSF_H(istate,jstate)-DegeneratePoint(ip).H(istate,jstate)
            end forall
            temp=dsyFrobeniusSquare(LSF_H,NStates)
            RMSDDegH=RMSDDegH+temp
            L=L+DegeneratePoint(ip).weight*(Ltemp+LSF_EnergyScaleSquare*temp)
            LSF_dcH=sy3UnitaryTransformation(dcHd_ByKnownf(LSF_f),LSF_phi,NExpansionCoefficients,NStates)
            do istate=1,NStates
                LSF_y( indicerow:indicerow+NStates-istate)=LSF_EnergyScaleSquare*DegeneratePoint(ip).H(istate:NStates,istate)
                LSF_MT(indicerow:indicerow+NStates-istate,:)=LSF_EnergyScaleSquare*transpose(LSF_dcH(:,istate:NStates,istate))
                indicerow=indicerow+NStates-istate+1
            end do
            !▽H (▽_c phi is neglected)
            LSF_dcdH=sy4UnitaryTransformation(dcdHd_ByKnownfdT(transpose(LSF_fd)),LSF_phi,NExpansionCoefficients,InternalDimension,NStates)
            do istate=1,NStates
                do jstate=istate,NStates
                    LSF_y( indicerow:indicerow+InternalDimension-1)=DegeneratePoint(ip).dH(:,jstate,istate)
                    LSF_MT(indicerow:indicerow+InternalDimension-1,:)=transpose(LSF_dcdH(:,:,jstate,istate))
                    indicerow=indicerow+InternalDimension
                end do
            end do
        end do
        RMSDDegH=Sqrt(RMSDDegH/(NStates*NStates)/NDegeneratePoints)
        RMSDDegdH=Sqrt(RMSDDegdH/(InternalDimension*NStates*NStates)/NDegeneratePoints)
        Ltemp=0d0
        do ip=1,NArtifactPoints!Unreliable data points, energy only
            call AdiabaticEnergy_State_f(ArtifactPoint(ip).geom,LSF_energy,LSF_phi,LSF_f)
            LSF_energy=LSF_energy-ArtifactPoint(ip).energy
            Ltemp=Ltemp+ArtifactPoint(ip).weight*dot_product(LSF_energy,LSF_energy)
            LSF_dcH=sy3UnitaryTransformation(dcHd_ByKnownf(LSF_f),LSF_phi,NExpansionCoefficients,NStates)
            LSF_y(indicerow:indicerow+NStates-1)=LSF_EnergyScaleSquare*ArtifactPoint(ip).energy
            forall(istate=1:NStates)
                LSF_MT(indicerow+istate-1,:)=LSF_EnergyScaleSquare*LSF_dcH(:,istate,istate)
            end forall
            indicerow=indicerow+NStates
        end do
        L=L+LSF_EnergyScaleSquare*Ltemp
        !Done construction, put them into A and b
        LSF_M=transpose(LSF_MT)
        forall(ip=1:LSF_NData)
            LSF_M(:,ip)=LSF_M(:,ip)*LSF_W(ip)
        end forall
        A=matmul(LSF_M,LSF_MT)
        forall(ip=1:NExpansionCoefficients)!Regularization
            A(ip,ip)=A(ip,ip)+LSF_Regularization
        end forall
        b=matmul(LSF_M,LSF_y)
    end subroutine LSFMatrices_L_RMSD
!---------------------------- End ----------------------------

!---------------- Trust region and dependency ----------------
    !In this section, to save memory, off-diagonals are treated as twice weighed
    !                 to save CPU time, weight -> Sqrt(weight)

    !Fit Hd by trust region method
    subroutine TrustRegion(c)
        real*8,dimension(NExpansionCoefficients),intent(inout)::c
        integer::i
        real*8::L,RMSDenergy,RMSDdH,RMSDDegH,RMSDDegdH
        !Initialize
            MaxMKLTrustRegionIteration=LSF_MaxLocalMinimizerIteration
            !weight -> Sqrt(weight)
                forall(i=1:NPoints)
                    point(i).weight=Sqrt(point(i).weight)
                end forall
                forall(i=1:NDegeneratePoints)
                    DegeneratePoint(i).weight=Sqrt(DegeneratePoint(i).weight)
                end forall
                forall(i=1:NArtifactPoints)
                    ArtifactPoint(i).weight=Sqrt(ArtifactPoint(i).weight)
                end forall
            !Allocate global local minimizer work space
                if(allocated(LSF_dHd)) deallocate(LSF_dHd)
                allocate(LSF_dHd(InternalDimension,NStates,NStates))
                if(allocated(LSF_dcphi)) deallocate(LSF_dcphi)
                allocate(LSF_dcphi(NExpansionCoefficients,NStates,NStates))
                if(allocated(LSF_dcHrep)) deallocate(LSF_dcHrep)
                allocate(LSF_dcHrep(NExpansionCoefficients,NStates,NStates))
                if(allocated(LSF_dcdHrep)) deallocate(LSF_dcdHrep)
                allocate(LSF_dcdHrep(NExpansionCoefficients,InternalDimension,NStates,NStates))
                !Trust region
                    if(allocated(LSF_spResidue)) deallocate(LSF_spResidue)
                    allocate(LSF_spResidue(DataPerPoint))
                    if(allocated(LSF_sdegpResidue)) deallocate(LSF_sdegpResidue)
                    allocate(LSF_sdegpResidue(DataPerDegeneratePoint))
                    if(allocated(LSF_spJacobian)) deallocate(LSF_spJacobian)
                    allocate(LSF_spJacobian(DataPerPoint,NExpansionCoefficients))
                    if(allocated(LSF_sdegpJacobian)) deallocate(LSF_sdegpJacobian)
                    allocate(LSF_sdegpJacobian(DataPerDegeneratePoint,NExpansionCoefficients))
        !Solve
        call showtime()
        write(*,'(1x,A40)')'Search for local minimum by trust region'
        call My_dtrnlsp(Residue,Jacobian,c,LSF_NData+NExpansionCoefficients,NExpansionCoefficients)
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
                deallocate(LSF_dHd)
                deallocate(LSF_dcphi)
                deallocate(LSF_dcHrep)
                deallocate(LSF_dcdHrep)
                !Trust region
                    deallocate(LSF_spResidue)
                    deallocate(LSF_sdegpResidue)
                    deallocate(LSF_spJacobian)
                    deallocate(LSF_sdegpJacobian)
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
    end subroutine TrustRegion

    subroutine Residue(r,c,M,N)
        integer,intent(in)::M,N
        real*8,dimension(M),intent(out)::r
        real*8,dimension(N),intent(in)::c
        integer::ip,istate,jstate,indicerow,indicesp
        real*8::Ltemp
        real*8,dimension(NStates,NStates)::phi
        !Initialize
            r(LSF_NData+1:M)=LSF_SqrtRegularization*c!Regularization
            call c2HdEC(c,HdEC,NExpansionCoefficients)
        indicerow=1!Start from 1st row
        do ip=1,NPoints!Regular data points
            call AdiabaticEnergy_dH(point(ip).geom,LSF_energy,LSF_dH)!Adiabatic representation
            call dFixdHPhase(LSF_dH,point(ip).dH,Ltemp,InternalDimension,NStates)!Fix off-diagonals phase
            LSF_spResidue(1:NStates)=LSF_EnergyScale*(LSF_energy-point(ip).energy)!Energy
            indicesp=NStates+1
            do istate=1,NStates!▽H
                LSF_spResidue(indicesp:indicesp+InternalDimension-1)=LSF_dH(:,istate,istate)-point(ip).dH(:,istate,istate)!Diagonal
                indicesp=indicesp+InternalDimension
                do jstate=istate+1,NStates!Off-diagonal
                    LSF_spResidue(indicesp:indicesp+InternalDimension-1)=Sqrt2*(LSF_dH(:,jstate,istate)-point(ip).dH(:,jstate,istate))
                    indicesp=indicesp+InternalDimension
                end do
            end do
            r(indicerow:indicerow+DataPerPoint-1)=point(ip).weight*LSF_spResidue
            indicerow=indicerow+DataPerPoint
        end do
        do ip=1,NDegeneratePoints!Almost degenerate data points
            call NondegenerateH_dH(DegeneratePoint(ip).geom,LSF_H,LSF_dH)!Nondegenerate representation
            call dFixHPhaseBydH(LSF_H,LSF_dH,DegeneratePoint(ip).dH,Ltemp,InternalDimension,NStates)!Fix off-diagonals phase
            indicesp=1
            do istate=1,NStates!H
                LSF_sdegpResidue(indicesp)=LSF_EnergyScale*(LSF_H(istate,istate)-DegeneratePoint(ip).H(istate,istate))!Diagonal
                LSF_sdegpResidue(indicesp+1:indicesp+NStates-istate)=Sqrt2*LSF_EnergyScale*(LSF_H(istate+1:NStates,istate)-DegeneratePoint(ip).H(istate+1:NStates,istate))
                indicesp=indicesp+NStates-istate+1
            end do
            do istate=1,NStates!▽H
                LSF_sdegpResidue(indicesp:indicesp+InternalDimension-1)=LSF_dH(:,istate,istate)-DegeneratePoint(ip).dH(:,istate,istate)!Diagonal
                indicesp=indicesp+InternalDimension
                do jstate=istate+1,NStates!Off-diagonal
                    LSF_sdegpResidue(indicesp:indicesp+InternalDimension-1)=Sqrt2*(LSF_dH(:,jstate,istate)-DegeneratePoint(ip).dH(:,jstate,istate))
                    indicesp=indicesp+InternalDimension
                end do
            end do
            r(indicerow:indicerow+DataPerDegeneratePoint-1)=DegeneratePoint(ip).weight*LSF_sdegpResidue
            indicerow=indicerow+DataPerDegeneratePoint
        end do
        do ip=1,NArtifactPoints!Unreliable data points, energy only
            r(indicerow:indicerow+NStates-1)=ArtifactPoint(ip).weight*LSF_EnergyScale*(AdiabaticEnergy(ArtifactPoint(ip).geom)-ArtifactPoint(ip).energy)
            indicerow=indicerow+NStates
        end do
    end subroutine Residue

    subroutine Jacobian(Jacob,c,M,N)
        integer,intent(in)::M,N
        real*8,dimension(M,N),intent(out)::Jacob
        real*8,dimension(N),intent(in)::c
        integer::ip,istate,jstate,i,indicerow,indicesp
        real*8::Ltemp
        real*8,dimension(NStates,NStates)::phi
        real*8,dimension(InternalDimension)::gradienttemp
        !Initialize
            Jacob(LSF_NData+1:M,:)=0d0
            forall(ip=1:NExpansionCoefficients)!Regularization
                Jacob(LSF_NData+ip,ip)=LSF_SqrtRegularization
            end forall
            call c2HdEC(c,HdEC,NExpansionCoefficients)
        indicerow=1!Start from 1st row
        do ip=1,NPoints!Regular data points
            call AdiabaticEnergy_dH_State_f_fd(point(ip).geom,LSF_energy,LSF_dH,LSF_phi,LSF_f,LSF_fd)!Adiabatic representation
            call dAssignBasisPhaseBydH(LSF_phi,LSF_dH,point(ip).dH,Ltemp,InternalDimension,NStates)!Assign basis phase
            !▽_c phi
            LSF_dcH=sy3UnitaryTransformation(dcHd_ByKnownf(LSF_f),LSF_phi,NExpansionCoefficients,NStates)
            LSF_dcphi=deigvec_ByKnowneigval_dA(LSF_energy,LSF_dcH,NExpansionCoefficients,NStates)
            !Energy
            forall(istate=1:NStates)
                LSF_spJacobian(istate,:)=LSF_EnergyScale*LSF_dcH(:,istate,istate)
            end forall
            !▽H
            LSF_dcdHrep=asy3matdirectmulsy3(LSF_dcphi,LSF_dH,NExpansionCoefficients,InternalDimension,NStates)
            LSF_dcdHrep=sy4UnitaryTransformation(dcdHd_ByKnownfdT(transpose(LSF_fd)),LSF_phi,NExpansionCoefficients,InternalDimension,NStates)&
                -LSF_dcdHrep-transpose4(LSF_dcdHrep,NExpansionCoefficients,InternalDimension,NStates,NStates)
            indicesp=NStates+1
            do istate=1,NStates
                LSF_spJacobian(indicesp:indicesp+InternalDimension-1,:)=transpose(LSF_dcdHrep(:,:,istate,istate))!Diagonal
                indicesp=indicesp+InternalDimension
                do jstate=istate+1,NStates!Off-diagonal
                    LSF_spJacobian(indicesp:indicesp+InternalDimension-1,:)=Sqrt2*transpose(LSF_dcdHrep(:,:,jstate,istate))
                    indicesp=indicesp+InternalDimension
                end do
            end do
            Jacob(indicerow:indicerow+DataPerPoint-1,:)=point(ip).weight*LSF_spJacobian
            indicerow=indicerow+DataPerPoint
        end do
        do ip=1,NDegeneratePoints!Almost degenerate data points
            !In this loop, LSF_energy stores the eigenvalues of nondegenerate operator
            call NondegenerateH_dH_eigval_State_dHd_f_fd(DegeneratePoint(ip).geom,LSF_H,LSF_dH,LSF_energy,LSF_phi,LSF_dHd,LSF_f,LSF_fd)!Nondegenerate representation
            call dFixHPhase_AssignBasisPhaseBydH(LSF_H,LSF_phi,LSF_dH,DegeneratePoint(ip).dH,Ltemp,InternalDimension,NStates)!Assign basis phase   
            !▽_c phi
            LSF_dcphi=deigvec_ByKnowneigval_dA(LSF_energy,sy3UnitaryTransformation(dcAd_ByKnown(dcdHd_ByKnownfdT(transpose(LSF_fd)),LSF_dHd),LSF_phi,NExpansionCoefficients,NStates),NExpansionCoefficients,NStates)
            !H
            LSF_dcHrep=asy3matmulsy(LSF_dcphi,LSF_H,NExpansionCoefficients,NStates)
            LSF_dcHrep=sy3UnitaryTransformation(dcHd_ByKnownf(LSF_f),LSF_phi,NExpansionCoefficients,NStates)&
                -LSF_dcHrep-transpose3(LSF_dcHrep,NExpansionCoefficients,NStates,NStates)
            indicesp=1
            do istate=1,NStates
                LSF_sdegpJacobian(indicesp,:)=LSF_EnergyScale*LSF_dcHrep(:,istate,istate)!Diagonal
                forall(jstate=istate+1:NStates)!Off-diagonal
                    LSF_sdegpJacobian(indicesp+jstate-istate,:)=Sqrt2*LSF_EnergyScale*LSF_dcHrep(:,jstate,istate)
                end forall
                indicesp=indicesp+NStates-istate+1
            end do
            !▽H
            LSF_dcdHrep=asy3matdirectmulsy3(LSF_dcphi,LSF_dH,NExpansionCoefficients,InternalDimension,NStates)
            LSF_dcdHrep=sy4UnitaryTransformation(dcdHd_ByKnownfdT(transpose(LSF_fd)),LSF_phi,NExpansionCoefficients,InternalDimension,NStates)&
                -LSF_dcdHrep-transpose4(LSF_dcdHrep,NExpansionCoefficients,InternalDimension,NStates,NStates)
            do istate=1,NStates
                LSF_spJacobian(indicesp:indicesp+InternalDimension-1,:)=transpose(LSF_dcdHrep(:,:,istate,istate))!Diagonal
                indicesp=indicesp+InternalDimension
                do jstate=istate+1,NStates!Off-diagonal
                    LSF_spJacobian(indicesp:indicesp+InternalDimension-1,:)=Sqrt2*transpose(LSF_dcdHrep(:,:,jstate,istate))
                    indicesp=indicesp+InternalDimension
                end do
            end do
            Jacob(indicerow:indicerow+DataPerDegeneratePoint-1,:)=DegeneratePoint(ip).weight*LSF_sdegpJacobian
            indicerow=indicerow+DataPerDegeneratePoint
        end do
        do ip=1,NArtifactPoints!Unreliable data points, energy only
            call AdiabaticEnergy_State_f(ArtifactPoint(ip).geom,LSF_energy,LSF_phi,LSF_f)
            LSF_dcH=sy3UnitaryTransformation(dcHd_ByKnownf(LSF_f),LSF_phi,NExpansionCoefficients,NStates)
            forall(istate=1:NStates)
                Jacob(indicerow+istate-1,:)=ArtifactPoint(ip).weight*LSF_EnergyScale*LSF_dcH(:,istate,istate)
            end forall
            indicerow=indicerow+NStates
        end do
    end subroutine Jacobian
!---------------------------- End ----------------------------

!---------- LBFGS/conjugate gradient and dependency ----------
    !To save CPU time, in this section Lagrangian -> Lagrangian / 2

    !Fit Hd by LBFGS method
    subroutine LimitedMemoryBFGS(c)
        real*8,dimension(NExpansionCoefficients),intent(inout)::c
        real*8::L,RMSDenergy,RMSDdH,RMSDDegH,RMSDDegdH
        !Initialize
            MaxQuasiNewtonIteration=LSF_MaxLocalMinimizerIteration
            !Allocate global local minimizer work space
                if(allocated(LSF_dHd)) deallocate(LSF_dHd)
                allocate(LSF_dHd(InternalDimension,NStates,NStates))
                if(allocated(LSF_dcphi)) deallocate(LSF_dcphi)
                allocate(LSF_dcphi(NExpansionCoefficients,NStates,NStates))
                if(allocated(LSF_dcHrep)) deallocate(LSF_dcHrep)
                allocate(LSF_dcHrep(NExpansionCoefficients,NStates,NStates))
                if(allocated(LSF_dcdHrep)) deallocate(LSF_dcdHrep)
                allocate(LSF_dcdHrep(NExpansionCoefficients,InternalDimension,NStates,NStates))
        !Solve
        call showtime()
        select case(LSF_LBFGSSolver)
            case('LBFGS')
                write(*,'(1x,A33)')'Search for local minimum by LBFGS'
                call LBFGS(Lagrangian,LagrangianGradient,c,NExpansionCoefficients,LSF_LBFGSMemory)
            case('LBFGS_Strong')
                write(*,'(1x,A62)')'Search for local minimum by LBFGS under strong Wolfe condition'
                call LBFGS_Strong_fdwithf(Lagrangian,LagrangianGradient,Lagrangian_LagrangianGradient,c,NExpansionCoefficients,LSF_LBFGSMemory)
            case default!Throw a warning
                write(*,'(1x,A39,1x,A32)')'Program abort: unsupported LBFGS solver',LSF_LBFGSSolver
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
            deallocate(LSF_dHd)
            deallocate(LSF_dcphi)
            deallocate(LSF_dcHrep)
            deallocate(LSF_dcdHrep)
    end subroutine LimitedMemoryBFGS

    !Fit Hd by conjugate gradient method
    subroutine ConjugateGradient(c)
        real*8,dimension(NExpansionCoefficients),intent(inout)::c
        real*8::L,RMSDenergy,RMSDdH,RMSDDegH,RMSDDegdH
        !Initialize
            MaxConjugateGradientIteration=LSF_MaxLocalMinimizerIteration
            !Allocate global local minimizer work space
                if(allocated(LSF_dHd)) deallocate(LSF_dHd)
                allocate(LSF_dHd(InternalDimension,NStates,NStates))
                if(allocated(LSF_dcphi)) deallocate(LSF_dcphi)
                allocate(LSF_dcphi(NExpansionCoefficients,NStates,NStates))
                if(allocated(LSF_dcHrep)) deallocate(LSF_dcHrep)
                allocate(LSF_dcHrep(NExpansionCoefficients,NStates,NStates))
                if(allocated(LSF_dcdHrep)) deallocate(LSF_dcdHrep)
                allocate(LSF_dcdHrep(NExpansionCoefficients,InternalDimension,NStates,NStates))
        !Solve
        call showtime()
        select case(LSF_CGSolver)
            case('DY')
                write(*,'(1x,A54)')'Search for local minimum by Dai-Yun conjugate gradient'
                call DYConjugateGradient(Lagrangian,LagrangianGradient,c,NExpansionCoefficients)
            case('DY_Strong')
                write(*,'(1x,A83)')'Search for local minimum by Dai-Yun conjugate gradient under strong Wolfe condition'
                call DYConjugateGradient_Strong_fdwithf(Lagrangian,LagrangianGradient,Lagrangian_LagrangianGradient,c,NExpansionCoefficients)
            case('PR')
                write(*,'(1x,A61)')'Search for local minimum by Polak-Ribiere+ conjugate gradient'
                call PRConjugateGradient_fdwithf(Lagrangian,LagrangianGradient,Lagrangian_LagrangianGradient,c,NExpansionCoefficients)
            case default!Throw a warning
                write(*,'(1x,A52,1x,A32)')'Program abort: unsupported conjugate gradient solver',LSF_CGSolver
                stop
        end select
        !Output
        call L_RMSD(c,L,RMSDenergy,RMSDdH,RMSDDegH,RMSDDegdH)
        write(*,'(1x,A29)')'Result of conjugate gradient:'
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
            deallocate(LSF_dHd)
            deallocate(LSF_dcphi)
            deallocate(LSF_dcHrep)
            deallocate(LSF_dcdHrep)
    end subroutine ConjugateGradient
    
    !dim dimensional vector c
    !L harvests Lagrangian
    subroutine Lagrangian(L,c,dim)
        real*8,intent(out)::L
        integer,intent(in)::dim
        real*8,dimension(dim),intent(in)::c
        integer::ip,istate,jstate
        real*8::Ltemp
        !Initialize
            L=LSF_Regularization*dot_product(c,c)!Regularization
            call c2HdEC(c,HdEC,NExpansionCoefficients)
        do ip=1,NPoints!Regular data points
            call AdiabaticEnergy_dH(point(ip).geom,LSF_energy,LSF_dH)!Adiabatic representation
            call dFixdHPhase(LSF_dH,point(ip).dH,Ltemp,InternalDimension,NStates)!Fix off-diagonals phase
            LSF_energy=LSF_energy-point(ip).energy!Energy (▽H has done during fixing)
            L=L+point(ip).weight*(Ltemp+LSF_EnergyScaleSquare*dot_product(LSF_energy,LSF_energy))
        end do
        do ip=1,NDegeneratePoints!Almost degenerate data points
            call NondegenerateH_dH(DegeneratePoint(ip).geom,LSF_H,LSF_dH)!Nondegenerate representation
            call dFixHPhaseBydH(LSF_H,LSF_dH,DegeneratePoint(ip).dH,Ltemp,InternalDimension,NStates)!Fix off-diagonals phase
            forall(istate=1:NStates,jstate=1:NStates,istate>=jstate)!H (▽H has done during fixing)
                LSF_H(istate,jstate)=LSF_H(istate,jstate)-DegeneratePoint(ip).H(istate,jstate)
            end forall
            L=L+DegeneratePoint(ip).weight*(Ltemp+LSF_EnergyScaleSquare*dsyFrobeniusSquare(LSF_H,NStates))
        end do
        Ltemp=0d0
        do ip=1,NArtifactPoints!Unreliable data points, energy only
            LSF_energy=AdiabaticEnergy(ArtifactPoint(ip).geom)-ArtifactPoint(ip).energy
            Ltemp=Ltemp+ArtifactPoint(ip).weight*dot_product(LSF_energy,LSF_energy)
        end do
        L=(L+LSF_EnergyScaleSquare*Ltemp)/2d0
    end subroutine Lagrangian

    !dim dimensional vector Ld, c
    !Ld harvests the gradient of Lagrangian over c
    subroutine LagrangianGradient(Ld,c,dim)
        integer,intent(in)::dim
        real*8,dimension(dim),intent(out)::Ld
        real*8,dimension(dim),intent(in)::c
        integer::i,ip,istate,jstate
        real*8::Ltemp
        real*8,dimension(dim)::Ldtemp
        !Initialize
            Ld=LSF_Regularization*c!Regularization
            call c2HdEC(c,HdEC,NExpansionCoefficients)
        do ip=1,NPoints!Regular data points
            call AdiabaticEnergy_dH_State_f_fd(point(ip).geom,LSF_energy,LSF_dH,LSF_phi,LSF_f,LSF_fd)!Adiabatic representation
            call dAssignBasisPhaseBydH(LSF_phi,LSF_dH,point(ip).dH,Ltemp,InternalDimension,NStates)!Assign basis phase
            !▽_c phi
            LSF_dcH=sy3UnitaryTransformation(dcHd_ByKnownf(LSF_f),LSF_phi,NExpansionCoefficients,NStates)
            LSF_dcphi=deigvec_ByKnowneigval_dA(LSF_energy,LSF_dcH,NExpansionCoefficients,NStates)
            !Energy
            forall(istate=1:NStates)
                LSF_dcHrep(:,istate,istate)=LSF_dcH(:,istate,istate)
            end forall
            forall(istate=2:NStates,jstate=1:NStates-1,istate>jstate)
                LSF_dcHrep(:,istate,jstate)=0d0
            end forall
            LSF_H=diag(LSF_energy-point(ip).energy,NStates)
            !▽H
            LSF_dcdHrep=asy3matdirectmulsy3(LSF_dcphi,LSF_dH,NExpansionCoefficients,InternalDimension,NStates)
            LSF_dcdHrep=sy4UnitaryTransformation(dcdHd_ByKnownfdT(transpose(LSF_fd)),LSF_phi,NExpansionCoefficients,InternalDimension,NStates)&
                -LSF_dcdHrep-transpose4(LSF_dcdHrep,NExpansionCoefficients,InternalDimension,NStates,NStates)
            forall(istate=1:NStates,jstate=1:NStates,istate>=jstate)
                LSF_dH(:,istate,jstate)=LSF_dH(:,istate,jstate)-point(ip).dH(:,istate,jstate)
            end forall
            Ld=Ld+point(ip).weight*(&
                LSF_EnergyScaleSquare*Trace3(sy3matmulsy(LSF_dcHrep,LSF_H,NExpansionCoefficients,NStates),NExpansionCoefficients,NStates)&    
                +Trace3(sy4matdotmulsy3(LSF_dcdHrep,LSF_dH,NExpansionCoefficients,InternalDimension,NStates),NExpansionCoefficients,NStates))
        end do
        do ip=1,NDegeneratePoints!Almost degenerate data points
            !In this loop, LSF_energy stores the eigenvalues of nondegenerate operator
            call NondegenerateH_dH_eigval_State_dHd_f_fd(DegeneratePoint(ip).geom,LSF_H,LSF_dH,LSF_energy,LSF_phi,LSF_dHd,LSF_f,LSF_fd)!Nondegenerate representation
            call dFixHPhase_AssignBasisPhaseBydH(LSF_H,LSF_phi,LSF_dH,DegeneratePoint(ip).dH,Ltemp,InternalDimension,NStates)!Assign basis phase
            !▽_c phi
            LSF_dcphi=deigvec_ByKnowneigval_dA(LSF_energy,sy3UnitaryTransformation(dcAd_ByKnown(dcdHd_ByKnownfdT(transpose(LSF_fd)),LSF_dHd),LSF_phi,NExpansionCoefficients,NStates),NExpansionCoefficients,NStates)
            !H
            LSF_dcHrep=asy3matmulsy(LSF_dcphi,LSF_H,NExpansionCoefficients,NStates)
            LSF_dcHrep=sy3UnitaryTransformation(dcHd_ByKnownf(LSF_f),LSF_phi,NExpansionCoefficients,NStates)&
                -LSF_dcHrep-transpose3(LSF_dcHrep,NExpansionCoefficients,NStates,NStates)
            forall(istate=1:NStates,jstate=1:NStates,istate>=jstate)
                LSF_H(istate,jstate)=LSF_H(istate,jstate)-DegeneratePoint(ip).H(istate,jstate)
            end forall
            !▽H
            LSF_dcdHrep=asy3matdirectmulsy3(LSF_dcphi,LSF_dH,NExpansionCoefficients,InternalDimension,NStates)
            LSF_dcdHrep=sy4UnitaryTransformation(dcdHd_ByKnownfdT(transpose(LSF_fd)),LSF_phi,NExpansionCoefficients,InternalDimension,NStates)&
                -LSF_dcdHrep-transpose4(LSF_dcdHrep,NExpansionCoefficients,InternalDimension,NStates,NStates)
            forall(istate=1:NStates,jstate=1:NStates,istate>=jstate)
                LSF_dH(:,istate,jstate)=LSF_dH(:,istate,jstate)-DegeneratePoint(ip).dH(:,istate,jstate)
            end forall
            Ld=Ld+DegeneratePoint(ip).weight*(&
                LSF_EnergyScaleSquare*Trace3(sy3matmulsy(LSF_dcHrep,LSF_H,NExpansionCoefficients,NStates),NExpansionCoefficients,NStates)&
                +Trace3(sy4matdotmulsy3(LSF_dcdHrep,LSF_dH,NExpansionCoefficients,InternalDimension,NStates),NExpansionCoefficients,NStates))
        end do
        Ldtemp=0d0
        do ip=1,NArtifactPoints!Unreliable data points, energy only
            call AdiabaticEnergy_State_f(ArtifactPoint(ip).geom,LSF_energy,LSF_phi,LSF_f)
            LSF_dcH=sy3UnitaryTransformation(dcHd_ByKnownf(LSF_f),LSF_phi,NExpansionCoefficients,NStates)
            forall(istate=1:NStates)
                LSF_dcHrep(:,istate,istate)=LSF_dcH(:,istate,istate)
            end forall
            forall(istate=2:NStates,jstate=1:NStates-1,istate>jstate)
                LSF_dcHrep(:,istate,jstate)=0d0
            end forall
            LSF_H=diag(LSF_energy-ArtifactPoint(ip).energy,NStates)
            Ldtemp=Ldtemp+ArtifactPoint(ip).weight*Trace3(sy3matmulsy(LSF_dcHrep,LSF_H,NExpansionCoefficients,NStates),NExpansionCoefficients,NStates)
        end do
        Ld=Ld+LSF_EnergyScaleSquare*Ldtemp
    end subroutine LagrangianGradient

    !dim dimensional vector Ld, c
    !L harvests Lagrangian
    !Ld harvests the gradient of Lagrangian over c
    subroutine Lagrangian_LagrangianGradient(L,Ld,c,dim)
        integer,intent(in)::dim
        real*8,intent(out)::L
        real*8,dimension(dim),intent(out)::Ld
        real*8,dimension(dim),intent(in)::c
        integer::i,ip,istate,jstate
        real*8::Ltemp
        real*8,dimension(dim)::Ldtemp
        !Initialize
            L =LSF_Regularization*dot_product(c,c)!Regularization
            Ld=LSF_Regularization*c!Regularization
            call c2HdEC(c,HdEC,NExpansionCoefficients)
        do ip=1,NPoints!Regular data points
            call AdiabaticEnergy_dH_State_f_fd(point(ip).geom,LSF_energy,LSF_dH,LSF_phi,LSF_f,LSF_fd)!Adiabatic representation
            call dAssignBasisPhaseBydH(LSF_phi,LSF_dH,point(ip).dH,Ltemp,InternalDimension,NStates)!Assign basis phase
            !▽_c phi
            LSF_dcH=sy3UnitaryTransformation(dcHd_ByKnownf(LSF_f),LSF_phi,NExpansionCoefficients,NStates)
            LSF_dcphi=deigvec_ByKnowneigval_dA(LSF_energy,LSF_dcH,NExpansionCoefficients,NStates)
            !Energy
            forall(istate=1:NStates)
                LSF_dcHrep(:,istate,istate)=LSF_dcH(:,istate,istate)
            end forall
            forall(istate=2:NStates,jstate=1:NStates-1,istate>jstate)
                LSF_dcHrep(:,istate,jstate)=0d0
            end forall
            LSF_H=diag(LSF_energy-point(ip).energy,NStates)
            !▽H
            LSF_dcdHrep=asy3matdirectmulsy3(LSF_dcphi,LSF_dH,NExpansionCoefficients,InternalDimension,NStates)
            LSF_dcdHrep=sy4UnitaryTransformation(dcdHd_ByKnownfdT(transpose(LSF_fd)),LSF_phi,NExpansionCoefficients,InternalDimension,NStates)&
                -LSF_dcdHrep-transpose4(LSF_dcdHrep,NExpansionCoefficients,InternalDimension,NStates,NStates)
            forall(istate=1:NStates,jstate=1:NStates,istate>=jstate)
                LSF_dH(:,istate,jstate)=LSF_dH(:,istate,jstate)-point(ip).dH(:,istate,jstate)
            end forall
            L =L +point(ip).weight*(Ltemp+LSF_EnergyScaleSquare*dsyFrobeniusSquare(LSF_H,NStates))
            Ld=Ld+point(ip).weight*(&
                LSF_EnergyScaleSquare*Trace3(sy3matmulsy(LSF_dcHrep,LSF_H,NExpansionCoefficients,NStates),NExpansionCoefficients,NStates)&    
                +Trace3(sy4matdotmulsy3(LSF_dcdHrep,LSF_dH,NExpansionCoefficients,InternalDimension,NStates),NExpansionCoefficients,NStates))
        end do
        do ip=1,NDegeneratePoints!Almost degenerate data points
            !In this loop, LSF_energy stores the eigenvalues of nondegenerate operator
            call NondegenerateH_dH_eigval_State_dHd_f_fd(DegeneratePoint(ip).geom,LSF_H,LSF_dH,LSF_energy,LSF_phi,LSF_dHd,LSF_f,LSF_fd)!Nondegenerate representation
            call dFixHPhase_AssignBasisPhaseBydH(LSF_H,LSF_phi,LSF_dH,DegeneratePoint(ip).dH,Ltemp,InternalDimension,NStates)!Assign basis phase
            !▽_c phi
            LSF_dcphi=deigvec_ByKnowneigval_dA(LSF_energy,sy3UnitaryTransformation(dcAd_ByKnown(dcdHd_ByKnownfdT(transpose(LSF_fd)),LSF_dHd),LSF_phi,NExpansionCoefficients,NStates),NExpansionCoefficients,NStates)
            !H
            LSF_dcHrep=asy3matmulsy(LSF_dcphi,LSF_H,NExpansionCoefficients,NStates)
            LSF_dcHrep=sy3UnitaryTransformation(dcHd_ByKnownf(LSF_f),LSF_phi,NExpansionCoefficients,NStates)&
                -LSF_dcHrep-transpose3(LSF_dcHrep,NExpansionCoefficients,NStates,NStates)
            forall(istate=1:NStates,jstate=1:NStates,istate>=jstate)
                LSF_H(istate,jstate)=LSF_H(istate,jstate)-DegeneratePoint(ip).H(istate,jstate)
            end forall
            !▽H
            LSF_dcdHrep=asy3matdirectmulsy3(LSF_dcphi,LSF_dH,NExpansionCoefficients,InternalDimension,NStates)
            LSF_dcdHrep=sy4UnitaryTransformation(dcdHd_ByKnownfdT(transpose(LSF_fd)),LSF_phi,NExpansionCoefficients,InternalDimension,NStates)&
                -LSF_dcdHrep-transpose4(LSF_dcdHrep,NExpansionCoefficients,InternalDimension,NStates,NStates)
            forall(istate=1:NStates,jstate=1:NStates,istate>=jstate)
                LSF_dH(:,istate,jstate)=LSF_dH(:,istate,jstate)-DegeneratePoint(ip).dH(:,istate,jstate)
            end forall
            L =L +DegeneratePoint(ip).weight*(Ltemp+LSF_EnergyScaleSquare*dsyFrobeniusSquare(LSF_H,NStates))
            Ld=Ld+DegeneratePoint(ip).weight*(&
                LSF_EnergyScaleSquare*Trace3(sy3matmulsy(LSF_dcHrep,LSF_H,NExpansionCoefficients,NStates),NExpansionCoefficients,NStates)&
                +Trace3(sy4matdotmulsy3(LSF_dcdHrep,LSF_dH,NExpansionCoefficients,InternalDimension,NStates),NExpansionCoefficients,NStates))
        end do
        Ltemp=0d0
        Ldtemp=0d0
        do ip=1,NArtifactPoints!Unreliable data points, energy only
            call AdiabaticEnergy_State_f(ArtifactPoint(ip).geom,LSF_energy,LSF_phi,LSF_f)
            LSF_dcH=sy3UnitaryTransformation(dcHd_ByKnownf(LSF_f),LSF_phi,NExpansionCoefficients,NStates)
            forall(istate=1:NStates)
                LSF_dcHrep(:,istate,istate)=LSF_dcH(:,istate,istate)
            end forall
            forall(istate=2:NStates,jstate=1:NStates-1,istate>jstate)
                LSF_dcHrep(:,istate,jstate)=0d0
            end forall
            LSF_H=diag(LSF_energy-ArtifactPoint(ip).energy,NStates)
            Ltemp = Ltemp+ArtifactPoint(ip).weight*dsyFrobeniusSquare(LSF_H,NStates)
            Ldtemp=Ldtemp+ArtifactPoint(ip).weight*Trace3(sy3matmulsy(LSF_dcHrep,LSF_H,NExpansionCoefficients,NStates),NExpansionCoefficients,NStates)
        end do
        L =(L+LSF_EnergyScaleSquare* Ltemp)/2d0
        Ld=Ld+LSF_EnergyScaleSquare*Ldtemp
    end subroutine Lagrangian_LagrangianGradient
!---------------------------- End ----------------------------

!--------------------- Auxiliary routine ---------------------
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
            L=LSF_Regularization*dot_product(c,c)!Regularization
            call c2HdEC(c,HdEC,NExpansionCoefficients)
        RMSDenergy=0d0
        RMSDdH=0d0
        do ip=1,NPoints!Regular data points, compute RMSD
            call AdiabaticEnergy_dH(point(ip).geom,LSF_energy,LSF_dH)!Adiabatic representation
            call dFixdHPhase(LSF_dH,point(ip).dH,Ltemp,InternalDimension,NStates)!Fix off-diagonals phase
            RMSDdH=RMSDdH+Ltemp
            LSF_energy=LSF_energy-point(ip).energy!Energy (dH has done during fixing)
            temp=dot_product(LSF_energy,LSF_energy)
            RMSDenergy=RMSDenergy+temp
            L=L+point(ip).weight*(Ltemp+LSF_EnergyScaleSquare*temp)
        end do
        RMSDenergy=Sqrt(RMSDenergy/NStates/NPoints)
        RMSDdH=Sqrt(RMSDdH/(InternalDimension*NStates*NStates)/NPoints)
        RMSDDegH=0d0
        RMSDDegdH=0d0
        do ip=1,NDegeneratePoints!Almost degenerate data points, compute RMSDDeg
            call NondegenerateH_dH(DegeneratePoint(ip).geom,LSF_H,LSF_dH)!Nondegenerate representation
            call dFixHPhaseBydH(LSF_H,LSF_dH,DegeneratePoint(ip).dH,Ltemp,InternalDimension,NStates)!Fix off-diagonals phase
            RMSDDegdH=RMSDDegdH+Ltemp
            forall(istate=1:NStates,jstate=1:NStates,istate>=jstate)!H (dH has done during fixing)
                LSF_H(istate,jstate)=LSF_H(istate,jstate)-DegeneratePoint(ip).H(istate,jstate)
            end forall
            temp=dsyFrobeniusSquare(LSF_H,NStates)
            RMSDDegH=RMSDDegH+temp
            L=L+DegeneratePoint(ip).weight*(Ltemp+LSF_EnergyScaleSquare*temp)
        end do
        RMSDDegH=Sqrt(RMSDDegH/(NStates*NStates)/NDegeneratePoints)
        RMSDDegdH=Sqrt(RMSDDegdH/(InternalDimension*NStates*NStates)/NDegeneratePoints)
        Ltemp=0d0
        do ip=1,NArtifactPoints!Unreliable data points, energy only
            LSF_energy=AdiabaticEnergy(ArtifactPoint(ip).geom)-ArtifactPoint(ip).energy
            Ltemp=Ltemp+ArtifactPoint(ip).weight*dot_product(LSF_energy,LSF_energy)
        end do
        L=L+LSF_EnergyScaleSquare*Ltemp
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
!---------------------------- End ----------------------------

end module HdLeastSquareFit