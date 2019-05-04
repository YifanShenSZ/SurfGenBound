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
        logical::HdLSF_UseStrongWolfe=.false.!Whether use strong Wolfe condition instead of Wolfe condition
        !LBFGS:
            integer::HdLSF_LBFGSMemory=10!Memory usage control, [3,30] is recommended (must > 0)
        !ConjugateGradient:
            character*32::HdLSF_ConjugateGradientSolver='DY'!Available: DY (Dai-Yun), PR (Polak-Ribiere+)

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
subroutine InitializeHdLeastSquareFit()!Initialize HdLeastSquareFit module
    integer::ip,istate,jstate
    real*8::MaxEnergy,MaxGrad,temp
    !A data point provides NState adiabatic energies, InternalDimension x NState x NState ▽H
    DataPerPoint=NState+NState*(NState+1)/2*InternalDimension!Upper triangle is redundant
    !A degenerate data point provides NState order H, InternalDimension x NState x NState ▽H
    DataPerDegeneratePoint=NState*(NState+1)/2*(InternalDimension+1)!Upper triangle is redundant
    MaxEnergy=0d0
    MaxGrad=0d0
    do ip=1,NPoints
        temp=maxval(Abs(point(ip).energy))
        if(temp>MaxEnergy) MaxEnergy=temp
        do istate=1,NState
            do jstate=istate,NState
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
        HdLSF_NData=DataPerPoint*NPoints+DataPerDegeneratePoint*NDegeneratePoints+NState*NArtifactPoints
        !Initial value of c
        allocate(c(NHdExpansionCoefficients))
        call HdEC2c(Hd_HdEC,c,NHdExpansionCoefficients)
        !Allocate global work space
            if(allocated(HdLSF_energy)) deallocate(HdLSF_energy)
            allocate(HdLSF_energy(NState))
            if(allocated(HdLSF_H)) deallocate(HdLSF_H)
            allocate(HdLSF_H(NState,NState))
            if(allocated(HdLSF_phi)) deallocate(HdLSF_phi)
            allocate(HdLSF_phi(NState,NState))
            if(allocated(HdLSF_dH)) deallocate(HdLSF_dH)
            allocate(HdLSF_dH(InternalDimension,NState,NState))
            if(allocated(HdLSF_dcH)) deallocate(HdLSF_dcH)
            allocate(HdLSF_dcH(NHdExpansionCoefficients,NState,NState))
            if(allocated(HdLSF_dcdH)) deallocate(HdLSF_dcdH)
            allocate(HdLSF_dcdH(NHdExpansionCoefficients,InternalDimension,NState,NState))
            if(allocated(HdLSF_f)) deallocate(HdLSF_f)
            allocate(HdLSF_f(NHdExpansionBasis))
            if(allocated(HdLSF_fd)) deallocate(HdLSF_fd)
            allocate(HdLSF_fd(InternalDimension,NHdExpansionBasis))
    !Fit Hd
    select case(HdLSF_Solver)
        !Single solvers
        case('pseudolinear')
            call pseudolinear(c)
        case('TrustRegion')
            call TrustRegionInterface(c)
        case('LineSearch')
            call LineSearchInterface(c)
        !2-step solvers
        case('pseudolinear_TrustRegion')
            do istate=1,HdLSF_Max2StepIteration
                call pseudolinear(c)
                call TrustRegionInterface(c)
            end do
        case('pseudolinear_LineSearch')
            do istate=1,HdLSF_Max2StepIteration
                call pseudolinear(c)
                call LineSearchInterface(c)
            end do
        case default!Throw a warning
            write(*,*)'Program abort: unsupported least square fit solver '//trim(adjustl(HdLSF_Solver))
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

!Lagrangian and root mean square deviations are the evaluation standards for the fit
    !Input:  current c
    !Output: L harvests Lagrangian,
    !        RMSDenergy/dH harvests root mean square deviation of adiabatic energy/dH over point,
    !        RMSDDegH/dH harvests root mean square deviation of nondegenerate H/dH over DegeneratePoint
subroutine L_RMSD(c,L,RMSDenergy,RMSDdH,RMSDDegH,RMSDDegdH)
	real*8,dimension(NHdExpansionCoefficients),intent(in)::c
	real*8,intent(out)::L,RMSDenergy,RMSDdH,RMSDDegH,RMSDDegdH
	integer::ip,istate,jstate
	real*8::Ltemp,temp
	!Initialize
		L=HdLSF_Regularization*dot_product(c,c)!Regularization
		call c2HdEC(c,Hd_HdEC,NHdExpansionCoefficients)
	RMSDenergy=0d0
	RMSDdH=0d0
	do ip=1,NPoints!Regular data points, compute RMSD
		call AdiabaticEnergy_dH(point(ip).geom,HdLSF_energy,HdLSF_dH)!Adiabatic representation
		call dFixdHPhase(HdLSF_dH,point(ip).dH,Ltemp,InternalDimension,NState)!Fix off-diagonals phase
		RMSDdH=RMSDdH+Ltemp
		HdLSF_energy=HdLSF_energy-point(ip).energy!Energy (dH has done during fixing)
		temp=dot_product(HdLSF_energy,HdLSF_energy)
		RMSDenergy=RMSDenergy+temp
		L=L+point(ip).weight*(Ltemp+HdLSF_EnergyScaleSquare*temp)
	end do
	RMSDenergy=dSqrt(RMSDenergy/NState/NPoints)
	RMSDdH=dSqrt(RMSDdH/(InternalDimension*NState*NState)/NPoints)
	RMSDDegH=0d0
	RMSDDegdH=0d0
	do ip=1,NDegeneratePoints!Almost degenerate data points, compute RMSDDeg
		call NondegenerateH_dH(DegeneratePoint(ip).geom,HdLSF_H,HdLSF_dH)!Nondegenerate representation
		call dFixHPhaseBydH(HdLSF_H,HdLSF_dH,DegeneratePoint(ip).dH,Ltemp,InternalDimension,NState)!Fix off-diagonals phase
		RMSDDegdH=RMSDDegdH+Ltemp
		forall(istate=1:NState,jstate=1:NState,istate>=jstate)!H (dH has done during fixing)
			HdLSF_H(istate,jstate)=HdLSF_H(istate,jstate)-DegeneratePoint(ip).H(istate,jstate)
		end forall
		temp=dsyFrobeniusSquare(HdLSF_H,NState)
		RMSDDegH=RMSDDegH+temp
		L=L+DegeneratePoint(ip).weight*(Ltemp+HdLSF_EnergyScaleSquare*temp)
	end do
	RMSDDegH=dSqrt(RMSDDegH/(NState*NState)/NDegeneratePoints)
	RMSDDegdH=dSqrt(RMSDDegdH/(InternalDimension*NState*NState)/NDegeneratePoints)
	Ltemp=0d0
	do ip=1,NArtifactPoints!Unreliable data points, energy only
		HdLSF_energy=AdiabaticEnergy(ArtifactPoint(ip).geom)-ArtifactPoint(ip).energy
		Ltemp=Ltemp+ArtifactPoint(ip).weight*dot_product(HdLSF_energy,HdLSF_energy)
	end do
	L=L+HdLSF_EnergyScaleSquare*Ltemp
end subroutine L_RMSD

!--------------- Solvers ---------------
    subroutine pseudolinear(cmin)!Fit Hd by pseudolinear method
        real*8,dimension(NHdExpansionCoefficients),intent(inout)::cmin
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
                        HdLSF_W(indice:indice+NState-1)=point(ip).weight!Energy
                        indice=indice+NState
                        do i=1,NState!▽H
                            HdLSF_W(indice:indice+InternalDimension-1)=point(ip).weight!Diagonal
                            indice=indice+InternalDimension
                            HdLSF_W(indice:indice+InternalDimension*(NState-i)-1)=2d0*point(ip).weight!Off-diagonal
                            indice=indice+InternalDimension*(NState-i)
                        end do
                    end do
                    do ip=1,NDegeneratePoints
                        do i=1,NState!H
                            HdLSF_W(indice)=DegeneratePoint(ip).weight!Diagonal
                            HdLSF_W(indice+1:indice+NState-i)=2d0*DegeneratePoint(ip).weight!Off-diagonal
                            indice=indice+NState-i+1
                        end do
                        do i=1,NState!▽H
                            HdLSF_W(indice:indice+InternalDimension-1)=DegeneratePoint(ip).weight!Diagonal
                            indice=indice+InternalDimension
                            HdLSF_W(indice:indice+InternalDimension*(NState-i)-1)=2d0*DegeneratePoint(ip).weight!Off-diagonal
                            indice=indice+InternalDimension*(NState-i)
                        end do
                    end do
                    do ip=1,NArtifactPoints
                        HdLSF_W(indice:indice+NState-1)=ArtifactPoint(ip).weight!Energy only
                        indice=indice+NState
                    end do
                if(allocated(HdLSF_M)) deallocate(HdLSF_M)
                allocate(HdLSF_M(NHdExpansionCoefficients,HdLSF_NData))
                if(allocated(HdLSF_MT)) deallocate(HdLSF_MT)
                allocate(HdLSF_MT(HdLSF_NData,NHdExpansionCoefficients))
            !Initialize linear least square fit and Lagrangian minimum
                allocate(c(NHdExpansionCoefficients))
                c=cmin
                allocate(b(NHdExpansionCoefficients))
                b=c
                allocate(A(NHdExpansionCoefficients,NHdExpansionCoefficients))
                call LSFMatrices_L(A,b,Lmin)
                L=Lmin
                indice=0
        !Solve
        call showtime()
        write(*,'(1x,A43)')'Explore phase space by pseudolinear hopping'
        do i=1,HdLSF_MaxHopperIteration!Main loop
            Lold=L!Prepare
            call My_dposv(A,b,NHdExpansionCoefficients)
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
        call c2HdEC(cmin,Hd_HdEC,NHdExpansionCoefficients)
        call WriteHdExpansionCoefficients()
        contains
            !The general form of weighted linear least square fit with pseudoregularization is:
            !    A c = b, where A = M . W . M^T + tau, b = M . W . y, c & tau have been explained at header,
            !    y is the data vector, M^T . c is the fitting prediction of y, W is the weight
            !To save memory, off-diagonals are treated as twice weighed
            !Input:  b = current c
            !Output: A harvests A, b harvests b, L harvests Lagrangian
            subroutine LSFMatrices_L(A,b,L)
                real*8,dimension(NHdExpansionCoefficients,NHdExpansionCoefficients),intent(out)::A
                real*8,dimension(NHdExpansionCoefficients),intent(inout)::b
                real*8,intent(out)::L
                integer::ip,istate,jstate,indicerow
                real*8::Ltemp
                !Initialize
                    L=HdLSF_Regularization*dot_product(b,b)!Regularization
                    call c2HdEC(b,Hd_HdEC,NHdExpansionCoefficients)
                !Construct M^T and y, add least square fit penalty to Lagrangian
                indicerow=1!Start from 1st row
                do ip=1,NPoints!Regular data points
                    call AdiabaticEnergy_dH_State_f_fd(point(ip).geom,HdLSF_energy,HdLSF_dH,HdLSF_phi,HdLSF_f,HdLSF_fd)!Adiabatic representation
                    call dAssignBasisPhaseBydH(HdLSF_phi,HdLSF_dH,point(ip).dH,Ltemp,InternalDimension,NState)!Assign basis phase
                    !Energy
                    HdLSF_energy=HdLSF_energy-point(ip).energy
                    L=L+point(ip).weight*(Ltemp+HdLSF_EnergyScaleSquare*dot_product(HdLSF_energy,HdLSF_energy))
                    HdLSF_dcH=sy3UnitaryTransformation(dcHd_ByKnownf(HdLSF_f),HdLSF_phi,NHdExpansionCoefficients,NState)
                    HdLSF_y(indicerow:indicerow+NState-1)=HdLSF_EnergyScaleSquare*point(ip).energy
                    forall(istate=1:NState)
                        HdLSF_MT(indicerow+istate-1,:)=HdLSF_EnergyScaleSquare*HdLSF_dcH(:,istate,istate)
                    end forall
                    indicerow=indicerow+NState
                    !▽H (▽_c phi is neglected)
                    HdLSF_dcdH=sy4UnitaryTransformation(dcdHd_ByKnownfdT(transpose(HdLSF_fd)),HdLSF_phi,NHdExpansionCoefficients,InternalDimension,NState)
                    do istate=1,NState
                        do jstate=istate,NState
                            HdLSF_y( indicerow:indicerow+InternalDimension-1)=point(ip).dH(:,jstate,istate)
                            HdLSF_MT(indicerow:indicerow+InternalDimension-1,:)=transpose(HdLSF_dcdH(:,:,jstate,istate))
                            indicerow=indicerow+InternalDimension
                        end do
                    end do
                end do
                do ip=1,NDegeneratePoints!Almost degenerate data points
                    call NondegenerateH_dH_State_f_fd(DegeneratePoint(ip).geom,HdLSF_H,HdLSF_dH,HdLSF_phi,HdLSF_f,HdLSF_fd)!Nondegenerate representation
                    call dFixHPhase_AssignBasisPhaseBydH(HdLSF_H,HdLSF_phi,HdLSF_dH,DegeneratePoint(ip).dH,Ltemp,InternalDimension,NState)!Assign basis phase
                    !H
                    forall(istate=1:NState,jstate=1:NState,istate>=jstate)
                        HdLSF_H(istate,jstate)=HdLSF_H(istate,jstate)-DegeneratePoint(ip).H(istate,jstate)
                    end forall
                    L=L+DegeneratePoint(ip).weight*(Ltemp+HdLSF_EnergyScaleSquare*dsyFrobeniusSquare(HdLSF_H,NState))
                    HdLSF_dcH=sy3UnitaryTransformation(dcHd_ByKnownf(HdLSF_f),HdLSF_phi,NHdExpansionCoefficients,NState)
                    do istate=1,NState
                        HdLSF_y( indicerow:indicerow+NState-istate)=HdLSF_EnergyScaleSquare*DegeneratePoint(ip).H(istate:NState,istate)
                        HdLSF_MT(indicerow:indicerow+NState-istate,:)=HdLSF_EnergyScaleSquare*transpose(HdLSF_dcH(:,istate:NState,istate))
                        indicerow=indicerow+NState-istate+1
                    end do
                    !▽H (▽_c phi is neglected)
                    HdLSF_dcdH=sy4UnitaryTransformation(dcdHd_ByKnownfdT(transpose(HdLSF_fd)),HdLSF_phi,NHdExpansionCoefficients,InternalDimension,NState)
                    do istate=1,NState
                        do jstate=istate,NState
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
                    HdLSF_dcH=sy3UnitaryTransformation(dcHd_ByKnownf(HdLSF_f),HdLSF_phi,NHdExpansionCoefficients,NState)
                    HdLSF_y(indicerow:indicerow+NState-1)=HdLSF_EnergyScaleSquare*ArtifactPoint(ip).energy
                    forall(istate=1:NState)
                        HdLSF_MT(indicerow+istate-1,:)=HdLSF_EnergyScaleSquare*HdLSF_dcH(:,istate,istate)
                    end forall
                    indicerow=indicerow+NState
                end do
                L=L+HdLSF_EnergyScaleSquare*Ltemp
                !Done construction, put them into A and b
                HdLSF_M=transpose(HdLSF_MT)
                forall(ip=1:HdLSF_NData)
                    HdLSF_M(:,ip)=HdLSF_M(:,ip)*HdLSF_W(ip)
                end forall
                A=matmul(HdLSF_M,HdLSF_MT)
                forall(ip=1:NHdExpansionCoefficients)!Regularization
                    A(ip,ip)=A(ip,ip)+HdLSF_Regularization
                end forall
                b=matmul(HdLSF_M,HdLSF_y)
            end subroutine LSFMatrices_L
            !Additional output:
            !    RMSDenergy/dH harvests root mean square deviation of adiabatic energy/dH over point,
            !    RMSDDegH/dH harvests root mean square deviation of nondegenerate H/dH over DegeneratePoint
            subroutine LSFMatrices_L_RMSD(A,b,L,RMSDenergy,RMSDdH,RMSDDegH,RMSDDegdH)
                real*8,dimension(NHdExpansionCoefficients,NHdExpansionCoefficients),intent(out)::A
                real*8,dimension(NHdExpansionCoefficients),intent(inout)::b
                real*8,intent(out)::L,RMSDenergy,RMSDdH,RMSDDegH,RMSDDegdH
                integer::ip,istate,jstate,indicerow
                real*8::Ltemp,temp
                !Initialize
                    L=HdLSF_Regularization*dot_product(b,b)!Regularization
                    call c2HdEC(b,Hd_HdEC,NHdExpansionCoefficients)
                !Construct M^T and y, add least square fit penalty to Lagrangian
                indicerow=1!Start from 1st row
                RMSDenergy=0d0
                RMSDdH=0d0
                do ip=1,NPoints!Regular data points, compute RMSD
                    call AdiabaticEnergy_dH_State_f_fd(point(ip).geom,HdLSF_energy,HdLSF_dH,HdLSF_phi,HdLSF_f,HdLSF_fd)!Adiabatic representation
                    call dAssignBasisPhaseBydH(HdLSF_phi,HdLSF_dH,point(ip).dH,Ltemp,InternalDimension,NState)!Assign basis phase
                    RMSDdH=RMSDdH+Ltemp
                    !Energy
                    HdLSF_energy=HdLSF_energy-point(ip).energy
                    temp=dot_product(HdLSF_energy,HdLSF_energy)
                    RMSDenergy=RMSDenergy+temp
                    L=L+point(ip).weight*(Ltemp+HdLSF_EnergyScaleSquare*temp)
                    HdLSF_dcH=sy3UnitaryTransformation(dcHd_ByKnownf(HdLSF_f),HdLSF_phi,NHdExpansionCoefficients,NState)
                    HdLSF_y(indicerow:indicerow+NState-1)=HdLSF_EnergyScaleSquare*point(ip).energy
                    forall(istate=1:NState)
                        HdLSF_MT(indicerow+istate-1,:)=HdLSF_EnergyScaleSquare*HdLSF_dcH(:,istate,istate)
                    end forall
                    indicerow=indicerow+NState
                    !▽H (▽_c phi is neglected)
                    HdLSF_dcdH=sy4UnitaryTransformation(dcdHd_ByKnownfdT(transpose(HdLSF_fd)),HdLSF_phi,NHdExpansionCoefficients,InternalDimension,NState)
                    do istate=1,NState
                        do jstate=istate,NState
                            HdLSF_y( indicerow:indicerow+InternalDimension-1)=point(ip).dH(:,jstate,istate)
                            HdLSF_MT(indicerow:indicerow+InternalDimension-1,:)=transpose(HdLSF_dcdH(:,:,jstate,istate))
                            indicerow=indicerow+InternalDimension
                        end do
                    end do
                end do
                RMSDenergy=dSqrt(RMSDenergy/NState/NPoints)
                RMSDdH=dSqrt(RMSDdH/(InternalDimension*NState*NState)/NPoints)
                RMSDDegH=0d0
                RMSDDegdH=0d0
                do ip=1,NDegeneratePoints!Almost degenerate data points, compute RMSDDeg
                    call NondegenerateH_dH_State_f_fd(DegeneratePoint(ip).geom,HdLSF_H,HdLSF_dH,HdLSF_phi,HdLSF_f,HdLSF_fd)!Nondegenerate representation
                    call dFixHPhase_AssignBasisPhaseBydH(HdLSF_H,HdLSF_phi,HdLSF_dH,DegeneratePoint(ip).dH,Ltemp,InternalDimension,NState)!Assign basis phase
                    RMSDDegdH=RMSDDegdH+Ltemp
                    !H
                    forall(istate=1:NState,jstate=1:NState,istate>=jstate)
                        HdLSF_H(istate,jstate)=HdLSF_H(istate,jstate)-DegeneratePoint(ip).H(istate,jstate)
                    end forall
                    temp=dsyFrobeniusSquare(HdLSF_H,NState)
                    RMSDDegH=RMSDDegH+temp
                    L=L+DegeneratePoint(ip).weight*(Ltemp+HdLSF_EnergyScaleSquare*temp)
                    HdLSF_dcH=sy3UnitaryTransformation(dcHd_ByKnownf(HdLSF_f),HdLSF_phi,NHdExpansionCoefficients,NState)
                    do istate=1,NState
                        HdLSF_y( indicerow:indicerow+NState-istate)=HdLSF_EnergyScaleSquare*DegeneratePoint(ip).H(istate:NState,istate)
                        HdLSF_MT(indicerow:indicerow+NState-istate,:)=HdLSF_EnergyScaleSquare*transpose(HdLSF_dcH(:,istate:NState,istate))
                        indicerow=indicerow+NState-istate+1
                    end do
                    !▽H (▽_c phi is neglected)
                    HdLSF_dcdH=sy4UnitaryTransformation(dcdHd_ByKnownfdT(transpose(HdLSF_fd)),HdLSF_phi,NHdExpansionCoefficients,InternalDimension,NState)
                    do istate=1,NState
                        do jstate=istate,NState
                            HdLSF_y( indicerow:indicerow+InternalDimension-1)=DegeneratePoint(ip).dH(:,jstate,istate)
                            HdLSF_MT(indicerow:indicerow+InternalDimension-1,:)=transpose(HdLSF_dcdH(:,:,jstate,istate))
                            indicerow=indicerow+InternalDimension
                        end do
                    end do
                end do
                RMSDDegH=dSqrt(RMSDDegH/(NState*NState)/NDegeneratePoints)
                RMSDDegdH=dSqrt(RMSDDegdH/(InternalDimension*NState*NState)/NDegeneratePoints)
                Ltemp=0d0
                do ip=1,NArtifactPoints!Unreliable data points, energy only
                    call AdiabaticEnergy_State_f(ArtifactPoint(ip).geom,HdLSF_energy,HdLSF_phi,HdLSF_f)
                    HdLSF_energy=HdLSF_energy-ArtifactPoint(ip).energy
                    Ltemp=Ltemp+ArtifactPoint(ip).weight*dot_product(HdLSF_energy,HdLSF_energy)
                    HdLSF_dcH=sy3UnitaryTransformation(dcHd_ByKnownf(HdLSF_f),HdLSF_phi,NHdExpansionCoefficients,NState)
                    HdLSF_y(indicerow:indicerow+NState-1)=HdLSF_EnergyScaleSquare*ArtifactPoint(ip).energy
                    forall(istate=1:NState)
                        HdLSF_MT(indicerow+istate-1,:)=HdLSF_EnergyScaleSquare*HdLSF_dcH(:,istate,istate)
                    end forall
                    indicerow=indicerow+NState
                end do
                L=L+HdLSF_EnergyScaleSquare*Ltemp
                !Done construction, put them into A and b
                HdLSF_M=transpose(HdLSF_MT)
                forall(ip=1:HdLSF_NData)
                    HdLSF_M(:,ip)=HdLSF_M(:,ip)*HdLSF_W(ip)
                end forall
                A=matmul(HdLSF_M,HdLSF_MT)
                forall(ip=1:NHdExpansionCoefficients)!Regularization
                    A(ip,ip)=A(ip,ip)+HdLSF_Regularization
                end forall
                b=matmul(HdLSF_M,HdLSF_y)
            end subroutine LSFMatrices_L_RMSD
    end subroutine pseudolinear
    
    subroutine TrustRegionInterface(c)!Fit Hd by trust region method
        !To save CPU time, weight -> Sqrt(weight), and will recover on exit
        real*8,dimension(NHdExpansionCoefficients),intent(inout)::c
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
                allocate(HdLSF_dHd(InternalDimension,NState,NState))
                if(allocated(HdLSF_dcphi)) deallocate(HdLSF_dcphi)
                allocate(HdLSF_dcphi(NHdExpansionCoefficients,NState,NState))
                if(allocated(HdLSF_dcHrep)) deallocate(HdLSF_dcHrep)
                allocate(HdLSF_dcHrep(NHdExpansionCoefficients,NState,NState))
                if(allocated(HdLSF_dcdHrep)) deallocate(HdLSF_dcdHrep)
                allocate(HdLSF_dcdHrep(NHdExpansionCoefficients,InternalDimension,NState,NState))
                !Trust region
                    if(allocated(HdLSF_spResidue)) deallocate(HdLSF_spResidue)
                    allocate(HdLSF_spResidue(DataPerPoint))
                    if(allocated(HdLSF_sdegpResidue)) deallocate(HdLSF_sdegpResidue)
                    allocate(HdLSF_sdegpResidue(DataPerDegeneratePoint))
                    if(allocated(HdLSF_spJacobian)) deallocate(HdLSF_spJacobian)
                    allocate(HdLSF_spJacobian(DataPerPoint,NHdExpansionCoefficients))
                    if(allocated(HdLSF_sdegpJacobian)) deallocate(HdLSF_sdegpJacobian)
                    allocate(HdLSF_sdegpJacobian(DataPerDegeneratePoint,NHdExpansionCoefficients))
        !Solve
        call showtime()
        write(*,'(1x,A40)')'Search for local minimum by trust region'
        call TrustRegion(Residue,c,HdLSF_NData+NHdExpansionCoefficients,NHdExpansionCoefficients,Jacobian=Jacobian,&
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
        call c2HdEC(c,Hd_HdEC,NHdExpansionCoefficients)
        call WriteHdExpansionCoefficients()
        contains
            !To save memory, off-diagonals are treated as twice weighed
            subroutine Residue(r,c,M,N)
                integer,intent(in)::M,N
                real*8,dimension(M),intent(out)::r
                real*8,dimension(N),intent(in)::c
                integer::ip,istate,jstate,indicerow,indicesp
                real*8::Ltemp
                real*8,dimension(NState,NState)::phi
                !Initialize
                    r(HdLSF_NData+1:M)=HdLSF_SqrtRegularization*c!Regularization
                    call c2HdEC(c,Hd_HdEC,NHdExpansionCoefficients)
                indicerow=1!Start from 1st row
                do ip=1,NPoints!Regular data points
                    call AdiabaticEnergy_dH(point(ip).geom,HdLSF_energy,HdLSF_dH)!Adiabatic representation
                    call dFixdHPhase(HdLSF_dH,point(ip).dH,Ltemp,InternalDimension,NState)!Fix off-diagonals phase
                    HdLSF_spResidue(1:NState)=HdLSF_EnergyScale*(HdLSF_energy-point(ip).energy)!Energy
                    indicesp=NState+1
                    do istate=1,NState!▽H
                        HdLSF_spResidue(indicesp:indicesp+InternalDimension-1)=HdLSF_dH(:,istate,istate)-point(ip).dH(:,istate,istate)!Diagonal
                        indicesp=indicesp+InternalDimension
                        do jstate=istate+1,NState!Off-diagonal
                            HdLSF_spResidue(indicesp:indicesp+InternalDimension-1)=Sqrt2*(HdLSF_dH(:,jstate,istate)-point(ip).dH(:,jstate,istate))
                            indicesp=indicesp+InternalDimension
                        end do
                    end do
                    r(indicerow:indicerow+DataPerPoint-1)=point(ip).weight*HdLSF_spResidue
                    indicerow=indicerow+DataPerPoint
                end do
                do ip=1,NDegeneratePoints!Almost degenerate data points
                    call NondegenerateH_dH(DegeneratePoint(ip).geom,HdLSF_H,HdLSF_dH)!Nondegenerate representation
                    call dFixHPhaseBydH(HdLSF_H,HdLSF_dH,DegeneratePoint(ip).dH,Ltemp,InternalDimension,NState)!Fix off-diagonals phase
                    indicesp=1
                    do istate=1,NState!H
                        HdLSF_sdegpResidue(indicesp)=HdLSF_EnergyScale*(HdLSF_H(istate,istate)-DegeneratePoint(ip).H(istate,istate))!Diagonal
                        HdLSF_sdegpResidue(indicesp+1:indicesp+NState-istate)=Sqrt2*HdLSF_EnergyScale*(HdLSF_H(istate+1:NState,istate)-DegeneratePoint(ip).H(istate+1:NState,istate))
                        indicesp=indicesp+NState-istate+1
                    end do
                    do istate=1,NState!▽H
                        HdLSF_sdegpResidue(indicesp:indicesp+InternalDimension-1)=HdLSF_dH(:,istate,istate)-DegeneratePoint(ip).dH(:,istate,istate)!Diagonal
                        indicesp=indicesp+InternalDimension
                        do jstate=istate+1,NState!Off-diagonal
                            HdLSF_sdegpResidue(indicesp:indicesp+InternalDimension-1)=Sqrt2*(HdLSF_dH(:,jstate,istate)-DegeneratePoint(ip).dH(:,jstate,istate))
                            indicesp=indicesp+InternalDimension
                        end do
                    end do
                    r(indicerow:indicerow+DataPerDegeneratePoint-1)=DegeneratePoint(ip).weight*HdLSF_sdegpResidue
                    indicerow=indicerow+DataPerDegeneratePoint
                end do
                do ip=1,NArtifactPoints!Unreliable data points, energy only
                    r(indicerow:indicerow+NState-1)=ArtifactPoint(ip).weight*HdLSF_EnergyScale*(AdiabaticEnergy(ArtifactPoint(ip).geom)-ArtifactPoint(ip).energy)
                    indicerow=indicerow+NState
                end do
            end subroutine Residue
            integer function Jacobian(Jacob,c,M,N)
                integer,intent(in)::M,N
                real*8,dimension(M,N),intent(out)::Jacob
                real*8,dimension(N),intent(in)::c
                integer::ip,istate,jstate,i,indicerow,indicesp
                real*8::Ltemp
                real*8,dimension(NState,NState)::phi
                real*8,dimension(InternalDimension)::gradienttemp
                !Initialize
                    Jacob(HdLSF_NData+1:M,:)=0d0
                    forall(ip=1:NHdExpansionCoefficients)!Regularization
                        Jacob(HdLSF_NData+ip,ip)=HdLSF_SqrtRegularization
                    end forall
                    call c2HdEC(c,Hd_HdEC,NHdExpansionCoefficients)
                indicerow=1!Start from 1st row
                do ip=1,NPoints!Regular data points
                    call AdiabaticEnergy_dH_State_f_fd(point(ip).geom,HdLSF_energy,HdLSF_dH,HdLSF_phi,HdLSF_f,HdLSF_fd)!Adiabatic representation
                    call dAssignBasisPhaseBydH(HdLSF_phi,HdLSF_dH,point(ip).dH,Ltemp,InternalDimension,NState)!Assign basis phase
                    !▽_c phi
                    HdLSF_dcH=sy3UnitaryTransformation(dcHd_ByKnownf(HdLSF_f),HdLSF_phi,NHdExpansionCoefficients,NState)
                    HdLSF_dcphi=deigvec_ByKnowneigval_dA(HdLSF_energy,HdLSF_dcH,NHdExpansionCoefficients,NState)
                    !Energy
                    forall(istate=1:NState)
                        HdLSF_spJacobian(istate,:)=HdLSF_EnergyScale*HdLSF_dcH(:,istate,istate)
                    end forall
                    !▽H
                    HdLSF_dcdHrep=asy3matdirectmulsy3(HdLSF_dcphi,HdLSF_dH,NHdExpansionCoefficients,InternalDimension,NState)
                    HdLSF_dcdHrep=sy4UnitaryTransformation(dcdHd_ByKnownfdT(transpose(HdLSF_fd)),HdLSF_phi,NHdExpansionCoefficients,InternalDimension,NState)&
                        -HdLSF_dcdHrep-transpose4(HdLSF_dcdHrep,NHdExpansionCoefficients,InternalDimension,NState,NState)
                    indicesp=NState+1
                    do istate=1,NState
                        HdLSF_spJacobian(indicesp:indicesp+InternalDimension-1,:)=transpose(HdLSF_dcdHrep(:,:,istate,istate))!Diagonal
                        indicesp=indicesp+InternalDimension
                        do jstate=istate+1,NState!Off-diagonal
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
                    call dFixHPhase_AssignBasisPhaseBydH(HdLSF_H,HdLSF_phi,HdLSF_dH,DegeneratePoint(ip).dH,Ltemp,InternalDimension,NState)!Assign basis phase   
                    !▽_c phi
                    HdLSF_dcphi=deigvec_ByKnowneigval_dA(HdLSF_energy,sy3UnitaryTransformation(dcAd_ByKnown(dcdHd_ByKnownfdT(transpose(HdLSF_fd)),HdLSF_dHd),HdLSF_phi,NHdExpansionCoefficients,NState),NHdExpansionCoefficients,NState)
                    !H
                    HdLSF_dcHrep=asy3matmulsy(HdLSF_dcphi,HdLSF_H,NHdExpansionCoefficients,NState)
                    HdLSF_dcHrep=sy3UnitaryTransformation(dcHd_ByKnownf(HdLSF_f),HdLSF_phi,NHdExpansionCoefficients,NState)&
                        -HdLSF_dcHrep-transpose3(HdLSF_dcHrep,NHdExpansionCoefficients,NState,NState)
                    indicesp=1
                    do istate=1,NState
                        HdLSF_sdegpJacobian(indicesp,:)=HdLSF_EnergyScale*HdLSF_dcHrep(:,istate,istate)!Diagonal
                        forall(jstate=istate+1:NState)!Off-diagonal
                            HdLSF_sdegpJacobian(indicesp+jstate-istate,:)=Sqrt2*HdLSF_EnergyScale*HdLSF_dcHrep(:,jstate,istate)
                        end forall
                        indicesp=indicesp+NState-istate+1
                    end do
                    !▽H
                    HdLSF_dcdHrep=asy3matdirectmulsy3(HdLSF_dcphi,HdLSF_dH,NHdExpansionCoefficients,InternalDimension,NState)
                    HdLSF_dcdHrep=sy4UnitaryTransformation(dcdHd_ByKnownfdT(transpose(HdLSF_fd)),HdLSF_phi,NHdExpansionCoefficients,InternalDimension,NState)&
                        -HdLSF_dcdHrep-transpose4(HdLSF_dcdHrep,NHdExpansionCoefficients,InternalDimension,NState,NState)
                    do istate=1,NState
                        HdLSF_spJacobian(indicesp:indicesp+InternalDimension-1,:)=transpose(HdLSF_dcdHrep(:,:,istate,istate))!Diagonal
                        indicesp=indicesp+InternalDimension
                        do jstate=istate+1,NState!Off-diagonal
                            HdLSF_spJacobian(indicesp:indicesp+InternalDimension-1,:)=Sqrt2*transpose(HdLSF_dcdHrep(:,:,jstate,istate))
                            indicesp=indicesp+InternalDimension
                        end do
                    end do
                    Jacob(indicerow:indicerow+DataPerDegeneratePoint-1,:)=DegeneratePoint(ip).weight*HdLSF_sdegpJacobian
                    indicerow=indicerow+DataPerDegeneratePoint
                end do
                do ip=1,NArtifactPoints!Unreliable data points, energy only
                    call AdiabaticEnergy_State_f(ArtifactPoint(ip).geom,HdLSF_energy,HdLSF_phi,HdLSF_f)
                    HdLSF_dcH=sy3UnitaryTransformation(dcHd_ByKnownf(HdLSF_f),HdLSF_phi,NHdExpansionCoefficients,NState)
                    forall(istate=1:NState)
                        Jacob(indicerow+istate-1,:)=ArtifactPoint(ip).weight*HdLSF_EnergyScale*HdLSF_dcH(:,istate,istate)
                    end forall
                    indicerow=indicerow+NState
                end do
                Jacobian=0!return 0
            end function Jacobian
    end subroutine TrustRegionInterface
    
    subroutine LineSearchInterface(c)!Fit Hd by line search method
        real*8,dimension(NHdExpansionCoefficients),intent(inout)::c
        real*8::L,RMSDenergy,RMSDdH,RMSDDegH,RMSDDegdH
        !Initialize
            !Allocate global local minimizer work space
            if(allocated(HdLSF_dHd)) deallocate(HdLSF_dHd)
            allocate(HdLSF_dHd(InternalDimension,NState,NState))
            if(allocated(HdLSF_dcphi)) deallocate(HdLSF_dcphi)
            allocate(HdLSF_dcphi(NHdExpansionCoefficients,NState,NState))
            if(allocated(HdLSF_dcHrep)) deallocate(HdLSF_dcHrep)
            allocate(HdLSF_dcHrep(NHdExpansionCoefficients,NState,NState))
            if(allocated(HdLSF_dcdHrep)) deallocate(HdLSF_dcdHrep)
            allocate(HdLSF_dcdHrep(NHdExpansionCoefficients,InternalDimension,NState,NState))
        !Solve
        call showtime()
        select case(HdLSF_LineSearcher)
            case('LBFGS')
                write(*,'(1x,A67)')'Search for local minimum by limited memory BFGS quasi-Newton method'
                call LBFGS(Lagrangian,LagrangianGradient,c,NHdExpansionCoefficients,f_fd=Lagrangian_LagrangianGradient,&
                    Memory=HdLSF_LBFGSMemory,Strong=HdLSF_UseStrongWolfe,MaxIteration=HdLSF_MaxLocalMinimizerIteration)
            case('ConjugateGradient')
                write(*,'(1x,A53)')'Search for local minimum by conjugate gradient method'
                call ConjugateGradient(Lagrangian,LagrangianGradient,c,NHdExpansionCoefficients,f_fd=Lagrangian_LagrangianGradient,&
                    Method=HdLSF_ConjugateGradientSolver,Strong=HdLSF_UseStrongWolfe,MaxIteration=HdLSF_MaxLocalMinimizerIteration)
            case default!Throw a warning
                write(*,*)'Program abort: unsupported line searcher '//trim(adjustl(HdLSF_LineSearcher))
                stop
        end select
        !Output
        call L_RMSD(c,L,RMSDenergy,RMSDdH,RMSDDegH,RMSDDegdH)
        write(*,'(1x,A22)')'Result of line search:'
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
        call c2HdEC(c,Hd_HdEC,NHdExpansionCoefficients)
        call WriteHdExpansionCoefficients()
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
                    call c2HdEC(c,Hd_HdEC,NHdExpansionCoefficients)
                do ip=1,NPoints!Regular data points
                    call AdiabaticEnergy_dH(point(ip).geom,HdLSF_energy,HdLSF_dH)!Adiabatic representation
                    call dFixdHPhase(HdLSF_dH,point(ip).dH,Ltemp,InternalDimension,NState)!Fix off-diagonals phase
                    HdLSF_energy=HdLSF_energy-point(ip).energy!Energy (▽H has done during fixing)
                    L=L+point(ip).weight*(Ltemp+HdLSF_EnergyScaleSquare*dot_product(HdLSF_energy,HdLSF_energy))
                end do
                do ip=1,NDegeneratePoints!Almost degenerate data points
                    call NondegenerateH_dH(DegeneratePoint(ip).geom,HdLSF_H,HdLSF_dH)!Nondegenerate representation
                    call dFixHPhaseBydH(HdLSF_H,HdLSF_dH,DegeneratePoint(ip).dH,Ltemp,InternalDimension,NState)!Fix off-diagonals phase
                    forall(istate=1:NState,jstate=1:NState,istate>=jstate)!H (▽H has done during fixing)
                        HdLSF_H(istate,jstate)=HdLSF_H(istate,jstate)-DegeneratePoint(ip).H(istate,jstate)
                    end forall
                    L=L+DegeneratePoint(ip).weight*(Ltemp+HdLSF_EnergyScaleSquare*dsyFrobeniusSquare(HdLSF_H,NState))
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
                    call c2HdEC(c,Hd_HdEC,NHdExpansionCoefficients)
                do ip=1,NPoints!Regular data points
                    call AdiabaticEnergy_dH_State_f_fd(point(ip).geom,HdLSF_energy,HdLSF_dH,HdLSF_phi,HdLSF_f,HdLSF_fd)!Adiabatic representation
                    call dAssignBasisPhaseBydH(HdLSF_phi,HdLSF_dH,point(ip).dH,Ltemp,InternalDimension,NState)!Assign basis phase
                    !▽_c phi
                    HdLSF_dcH=sy3UnitaryTransformation(dcHd_ByKnownf(HdLSF_f),HdLSF_phi,NHdExpansionCoefficients,NState)
                    HdLSF_dcphi=deigvec_ByKnowneigval_dA(HdLSF_energy,HdLSF_dcH,NHdExpansionCoefficients,NState)
                    !Energy
                    forall(istate=1:NState)
                        HdLSF_dcHrep(:,istate,istate)=HdLSF_dcH(:,istate,istate)
                    end forall
                    forall(istate=2:NState,jstate=1:NState-1,istate>jstate)
                        HdLSF_dcHrep(:,istate,jstate)=0d0
                    end forall
                    HdLSF_H=diag(HdLSF_energy-point(ip).energy,NState)
                    !▽H
                    HdLSF_dcdHrep=asy3matdirectmulsy3(HdLSF_dcphi,HdLSF_dH,NHdExpansionCoefficients,InternalDimension,NState)
                    HdLSF_dcdHrep=sy4UnitaryTransformation(dcdHd_ByKnownfdT(transpose(HdLSF_fd)),HdLSF_phi,NHdExpansionCoefficients,InternalDimension,NState)&
                        -HdLSF_dcdHrep-transpose4(HdLSF_dcdHrep,NHdExpansionCoefficients,InternalDimension,NState,NState)
                    forall(istate=1:NState,jstate=1:NState,istate>=jstate)
                        HdLSF_dH(:,istate,jstate)=HdLSF_dH(:,istate,jstate)-point(ip).dH(:,istate,jstate)
                    end forall
                    Ld=Ld+point(ip).weight*(&
                        HdLSF_EnergyScaleSquare*Trace3(sy3matmulsy(HdLSF_dcHrep,HdLSF_H,NHdExpansionCoefficients,NState),NHdExpansionCoefficients,NState)&    
                        +Trace3(sy4matdotmulsy3(HdLSF_dcdHrep,HdLSF_dH,NHdExpansionCoefficients,InternalDimension,NState),NHdExpansionCoefficients,NState))
                end do
                do ip=1,NDegeneratePoints!Almost degenerate data points
                    !In this loop, HdLSF_energy stores the eigenvalues of nondegenerate operator
                    call NondegenerateH_dH_eigval_State_dHd_f_fd(DegeneratePoint(ip).geom,HdLSF_H,HdLSF_dH,HdLSF_energy,HdLSF_phi,HdLSF_dHd,HdLSF_f,HdLSF_fd)!Nondegenerate representation
                    call dFixHPhase_AssignBasisPhaseBydH(HdLSF_H,HdLSF_phi,HdLSF_dH,DegeneratePoint(ip).dH,Ltemp,InternalDimension,NState)!Assign basis phase
                    !▽_c phi
                    HdLSF_dcphi=deigvec_ByKnowneigval_dA(HdLSF_energy,sy3UnitaryTransformation(dcAd_ByKnown(dcdHd_ByKnownfdT(transpose(HdLSF_fd)),HdLSF_dHd),HdLSF_phi,NHdExpansionCoefficients,NState),NHdExpansionCoefficients,NState)
                    !H
                    HdLSF_dcHrep=asy3matmulsy(HdLSF_dcphi,HdLSF_H,NHdExpansionCoefficients,NState)
                    HdLSF_dcHrep=sy3UnitaryTransformation(dcHd_ByKnownf(HdLSF_f),HdLSF_phi,NHdExpansionCoefficients,NState)&
                        -HdLSF_dcHrep-transpose3(HdLSF_dcHrep,NHdExpansionCoefficients,NState,NState)
                    forall(istate=1:NState,jstate=1:NState,istate>=jstate)
                        HdLSF_H(istate,jstate)=HdLSF_H(istate,jstate)-DegeneratePoint(ip).H(istate,jstate)
                    end forall
                    !▽H
                    HdLSF_dcdHrep=asy3matdirectmulsy3(HdLSF_dcphi,HdLSF_dH,NHdExpansionCoefficients,InternalDimension,NState)
                    HdLSF_dcdHrep=sy4UnitaryTransformation(dcdHd_ByKnownfdT(transpose(HdLSF_fd)),HdLSF_phi,NHdExpansionCoefficients,InternalDimension,NState)&
                        -HdLSF_dcdHrep-transpose4(HdLSF_dcdHrep,NHdExpansionCoefficients,InternalDimension,NState,NState)
                    forall(istate=1:NState,jstate=1:NState,istate>=jstate)
                        HdLSF_dH(:,istate,jstate)=HdLSF_dH(:,istate,jstate)-DegeneratePoint(ip).dH(:,istate,jstate)
                    end forall
                    Ld=Ld+DegeneratePoint(ip).weight*(&
                        HdLSF_EnergyScaleSquare*Trace3(sy3matmulsy(HdLSF_dcHrep,HdLSF_H,NHdExpansionCoefficients,NState),NHdExpansionCoefficients,NState)&
                        +Trace3(sy4matdotmulsy3(HdLSF_dcdHrep,HdLSF_dH,NHdExpansionCoefficients,InternalDimension,NState),NHdExpansionCoefficients,NState))
                end do
                Ldtemp=0d0
                do ip=1,NArtifactPoints!Unreliable data points, energy only
                    call AdiabaticEnergy_State_f(ArtifactPoint(ip).geom,HdLSF_energy,HdLSF_phi,HdLSF_f)
                    HdLSF_dcH=sy3UnitaryTransformation(dcHd_ByKnownf(HdLSF_f),HdLSF_phi,NHdExpansionCoefficients,NState)
                    forall(istate=1:NState)
                        HdLSF_dcHrep(:,istate,istate)=HdLSF_dcH(:,istate,istate)
                    end forall
                    forall(istate=2:NState,jstate=1:NState-1,istate>jstate)
                        HdLSF_dcHrep(:,istate,jstate)=0d0
                    end forall
                    HdLSF_H=diag(HdLSF_energy-ArtifactPoint(ip).energy,NState)
                    Ldtemp=Ldtemp+ArtifactPoint(ip).weight*Trace3(sy3matmulsy(HdLSF_dcHrep,HdLSF_H,NHdExpansionCoefficients,NState),NHdExpansionCoefficients,NState)
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
                    call c2HdEC(c,Hd_HdEC,NHdExpansionCoefficients)
                do ip=1,NPoints!Regular data points
                    call AdiabaticEnergy_dH_State_f_fd(point(ip).geom,HdLSF_energy,HdLSF_dH,HdLSF_phi,HdLSF_f,HdLSF_fd)!Adiabatic representation
                    call dAssignBasisPhaseBydH(HdLSF_phi,HdLSF_dH,point(ip).dH,Ltemp,InternalDimension,NState)!Assign basis phase
                    !▽_c phi
                    HdLSF_dcH=sy3UnitaryTransformation(dcHd_ByKnownf(HdLSF_f),HdLSF_phi,NHdExpansionCoefficients,NState)
                    HdLSF_dcphi=deigvec_ByKnowneigval_dA(HdLSF_energy,HdLSF_dcH,NHdExpansionCoefficients,NState)
                    !Energy
                    forall(istate=1:NState)
                        HdLSF_dcHrep(:,istate,istate)=HdLSF_dcH(:,istate,istate)
                    end forall
                    forall(istate=2:NState,jstate=1:NState-1,istate>jstate)
                        HdLSF_dcHrep(:,istate,jstate)=0d0
                    end forall
                    HdLSF_H=diag(HdLSF_energy-point(ip).energy,NState)
                    !▽H
                    HdLSF_dcdHrep=asy3matdirectmulsy3(HdLSF_dcphi,HdLSF_dH,NHdExpansionCoefficients,InternalDimension,NState)
                    HdLSF_dcdHrep=sy4UnitaryTransformation(dcdHd_ByKnownfdT(transpose(HdLSF_fd)),HdLSF_phi,NHdExpansionCoefficients,InternalDimension,NState)&
                        -HdLSF_dcdHrep-transpose4(HdLSF_dcdHrep,NHdExpansionCoefficients,InternalDimension,NState,NState)
                    forall(istate=1:NState,jstate=1:NState,istate>=jstate)
                        HdLSF_dH(:,istate,jstate)=HdLSF_dH(:,istate,jstate)-point(ip).dH(:,istate,jstate)
                    end forall
                    L =L +point(ip).weight*(Ltemp+HdLSF_EnergyScaleSquare*dsyFrobeniusSquare(HdLSF_H,NState))
                    Ld=Ld+point(ip).weight*(&
                        HdLSF_EnergyScaleSquare*Trace3(sy3matmulsy(HdLSF_dcHrep,HdLSF_H,NHdExpansionCoefficients,NState),NHdExpansionCoefficients,NState)&    
                        +Trace3(sy4matdotmulsy3(HdLSF_dcdHrep,HdLSF_dH,NHdExpansionCoefficients,InternalDimension,NState),NHdExpansionCoefficients,NState))
                end do
                do ip=1,NDegeneratePoints!Almost degenerate data points
                    !In this loop, HdLSF_energy stores the eigenvalues of nondegenerate operator
                    call NondegenerateH_dH_eigval_State_dHd_f_fd(DegeneratePoint(ip).geom,HdLSF_H,HdLSF_dH,HdLSF_energy,HdLSF_phi,HdLSF_dHd,HdLSF_f,HdLSF_fd)!Nondegenerate representation
                    call dFixHPhase_AssignBasisPhaseBydH(HdLSF_H,HdLSF_phi,HdLSF_dH,DegeneratePoint(ip).dH,Ltemp,InternalDimension,NState)!Assign basis phase
                    !▽_c phi
                    HdLSF_dcphi=deigvec_ByKnowneigval_dA(HdLSF_energy,sy3UnitaryTransformation(dcAd_ByKnown(dcdHd_ByKnownfdT(transpose(HdLSF_fd)),HdLSF_dHd),HdLSF_phi,NHdExpansionCoefficients,NState),NHdExpansionCoefficients,NState)
                    !H
                    HdLSF_dcHrep=asy3matmulsy(HdLSF_dcphi,HdLSF_H,NHdExpansionCoefficients,NState)
                    HdLSF_dcHrep=sy3UnitaryTransformation(dcHd_ByKnownf(HdLSF_f),HdLSF_phi,NHdExpansionCoefficients,NState)&
                        -HdLSF_dcHrep-transpose3(HdLSF_dcHrep,NHdExpansionCoefficients,NState,NState)
                    forall(istate=1:NState,jstate=1:NState,istate>=jstate)
                        HdLSF_H(istate,jstate)=HdLSF_H(istate,jstate)-DegeneratePoint(ip).H(istate,jstate)
                    end forall
                    !▽H
                    HdLSF_dcdHrep=asy3matdirectmulsy3(HdLSF_dcphi,HdLSF_dH,NHdExpansionCoefficients,InternalDimension,NState)
                    HdLSF_dcdHrep=sy4UnitaryTransformation(dcdHd_ByKnownfdT(transpose(HdLSF_fd)),HdLSF_phi,NHdExpansionCoefficients,InternalDimension,NState)&
                        -HdLSF_dcdHrep-transpose4(HdLSF_dcdHrep,NHdExpansionCoefficients,InternalDimension,NState,NState)
                    forall(istate=1:NState,jstate=1:NState,istate>=jstate)
                        HdLSF_dH(:,istate,jstate)=HdLSF_dH(:,istate,jstate)-DegeneratePoint(ip).dH(:,istate,jstate)
                    end forall
                    L =L +DegeneratePoint(ip).weight*(Ltemp+HdLSF_EnergyScaleSquare*dsyFrobeniusSquare(HdLSF_H,NState))
                    Ld=Ld+DegeneratePoint(ip).weight*(&
                        HdLSF_EnergyScaleSquare*Trace3(sy3matmulsy(HdLSF_dcHrep,HdLSF_H,NHdExpansionCoefficients,NState),NHdExpansionCoefficients,NState)&
                        +Trace3(sy4matdotmulsy3(HdLSF_dcdHrep,HdLSF_dH,NHdExpansionCoefficients,InternalDimension,NState),NHdExpansionCoefficients,NState))
                end do
                Ltemp=0d0
                Ldtemp=0d0
                do ip=1,NArtifactPoints!Unreliable data points, energy only
                    call AdiabaticEnergy_State_f(ArtifactPoint(ip).geom,HdLSF_energy,HdLSF_phi,HdLSF_f)
                    HdLSF_dcH=sy3UnitaryTransformation(dcHd_ByKnownf(HdLSF_f),HdLSF_phi,NHdExpansionCoefficients,NState)
                    forall(istate=1:NState)
                        HdLSF_dcHrep(:,istate,istate)=HdLSF_dcH(:,istate,istate)
                    end forall
                    forall(istate=2:NState,jstate=1:NState-1,istate>jstate)
                        HdLSF_dcHrep(:,istate,jstate)=0d0
                    end forall
                    HdLSF_H=diag(HdLSF_energy-ArtifactPoint(ip).energy,NState)
                    Ltemp = Ltemp+ArtifactPoint(ip).weight*dsyFrobeniusSquare(HdLSF_H,NState)
                    Ldtemp=Ldtemp+ArtifactPoint(ip).weight*Trace3(sy3matmulsy(HdLSF_dcHrep,HdLSF_H,NHdExpansionCoefficients,NState),NHdExpansionCoefficients,NState)
                end do
                L =(L+HdLSF_EnergyScaleSquare* Ltemp)/2d0
                Ld=Ld+HdLSF_EnergyScaleSquare*Ldtemp
                Lagrangian_LagrangianGradient=0!return 0
            end function Lagrangian_LagrangianGradient
    end subroutine LineSearchInterface
!----------------- End -----------------

end module HdLeastSquareFit