!Provide analyzation of the fitted potential energy surface:
!    Compute H & ▽H in diabatic & adiabatic & nondegenerate representation on specified geometries
!    Minimum search and vibration analysis of specified adiabatic state
!    Mex search and gh orthogonalization along conical intersection seam between specified adiabatic states
module Analyzation
    use Basic
    implicit none

!Analyzation module only variable
	!Analyzation input
        character*32::AnalyzationJobType
		integer::InterestingState
	!Geometry input
		integer::NGeoms
		real*8,allocatable,dimension(:,:)::InterestingGeom
		real*8,allocatable,dimension(:)::g_AbInitio,h_AbInitio
    
contains
subroutine Analyze()!Top level standard interface for other modules to call
	call ReadAnalyzeInput()
    select case(AnalyzationJobType)
        case('min')
            call MinimumSearch()
        case('mex')
            call MexSearch()
        case default!Throw a warning
            write(*,'(1x,A47,1x,A32)')'Program abort: unsupported analyzation job type',AnalyzationJobType
            stop
    end select
end subroutine Analyze
!BUG TO FIX: Where q is the internal coordinate (for this program it's internal coordinate difference)
subroutine ReadAnalyzeInput()!Read the input file for Analyzation: AnalyzeInput
	character*128::InterestingGeomFile
	integer::i
	real*8,dimension(3)::buffer
    open(unit=99,file='analyzation.in',status='old')
        read(99,*)
        read(99,*)
        read(99,*)
		read(99,*)AnalyzationJobType
		    write(*,*)'Analyzation job type: '//AnalyzationJobType
        read(99,*)
		read(99,*)InterestingState
		read(99,*)
        read(99,*)InterestingGeomFile
	close(99)
	open(unit=99,file=InterestingGeomFile,status='old')
	    NGeoms=0!Count number of geometries
		do
			read(99,*,iostat=i)buffer
			if(i/=0) exit
			NGeoms=NGeoms+1
		end do
		NGeoms=NGeoms/NAtoms
		rewind 99
		allocate(InterestingGeom(CartesianDimension,NGeoms))!Read geometries
		do i=1,NGeoms
			read(99,*)InterestingGeom(:,i)
		end do
	close(99)
end subroutine ReadAnalyzeInput

subroutine MinimumSearch()
    integer::i
    real*8,dimension(InternalDimension)::q,freq
    real*8,dimension(CartesianDimension)::r
    real*8,dimension(InternalDimension,CartesianDimension)::B
	real*8,dimension(InternalDimension,InternalDimension)::Hessian
	write(*,'(1x,A46,1x,I2)')'Search for minimum on potential energy surface',InterestingState
    q=InternalCoordinateq(InterestingGeom(:,1),InternalDImension,CartesianDimension)
    call BFGS(AdiabaticEnergyInterface,AdiabaticGradientInterface,q,InternalDImension,&
        fdd=AdiabaticHessianInterface,f_fd=AdiabaticEnergy_GradientInterface)
    r=CartesianCoordinater(q,CartesianDimension,InternalDImension,MoleculeDetail.mass,InterestingGeom(:,1))
    call WilsonBMatrixAndInternalCoordinateq(B,q,r,InternalDImension,CartesianDimension)
    i=AdiabaticHessianInterface(Hessian,q,InternalDimension)
    call WilsonGFMethod(freq,Hessian,InternalDimension,B,MoleculeDetail.mass,NAtoms)
	open(unit=99,file='MinimumCartesianGeometry.xyz',status='replace')
		write(99,*)NAtoms
		write(99,*)
        do i=1,NAtoms
            write(99,'(A2,3F20.15)')MoleculeDetail.ElementSymbol(i),r(3*i-2:3*i)
        end do
    close(99)
    open(unit=99,file='MinimumInternalGeometry.out',status='replace')
        write(99,*)q
    close(99)
    open(unit=99,file='MinimumVibrationalFrequency.txt',status='replace')
        write(99,'(A19)')'Mode'//char(9)//'Frequency/cm-1'
        do i=1,InternalDimension
            write(99,'(A4,1x,F14.8)')i,freq(i)/cm_1InAu
        end do
    close(99)
    open(unit=99,file='MinimumNormalMode.out',status='replace')
        write(99,*)Hessian
    close(99)
end subroutine MinimumSearch

subroutine MexSearch()
    integer::i
    real*8,dimension(InternalDimension)::q
	real*8,dimension(CartesianDimension)::r
	real*8,dimension(InternalDimension,NStates,NStates)::dH
	write(*,'(1x,A48,1x,I2,1x,A3,1x,I2)')'Search for mex between potential energy surfaces',InterestingState,'and',InterestingState+1
    q=InternalCoordinateq(InterestingGeom(:,1),InternalDImension,CartesianDimension)
    if(NStates==2) then!2 state case we can simply search for minimum of Hd diagonal subject to zero off-diagonal and degenerate diagonals
        call AugmentedLagrangian(f,fd,c,cd,q,InternalDImension,2,&
            fdd=fdd,cdd=cdd,Precision=1d-8)!This is Columbus7 energy precision
    else!In general case we have to search for minimum on potential energy surface of interest subject to degeneracy constaint
	    call AugmentedLagrangian(AdiabaticEnergyInterface,AdiabaticGradientInterface,AdiabaticGapInterface,AdiabaticGapGradientInterface,&
	        q,InternalDImension,1,Precision=1d-8,&!This is Columbus7 energy precision
            fdd=AdiabaticHessianInterface,cdd=AdiabaticGapHessianInterface,f_fd=AdiabaticEnergy_GradientInterface)
    end if
	r=CartesianCoordinater(q,CartesianDimension,InternalDImension,MoleculeDetail.mass,InterestingGeom(:,1))
	dH=AdiabaticdH(q)
	if(allocated(g_AbInitio).and.allocated(h_AbInitio)) then
		call ghOrthogonalization(dH(:,InterestingState,InterestingState),dH(:,InterestingState+1,InterestingState+1),dH(:,InterestingState+1,InterestingState),InternalDimension,&
		    gref=g_AbInitio,href=h_AbInitio)
	else
	    call ghOrthogonalization(dH(:,InterestingState,InterestingState),dH(:,InterestingState+1,InterestingState+1),dH(:,InterestingState+1,InterestingState),InternalDimension)
	end if
	open(unit=99,file='MexCartesianGeometry.xyz',status='replace')
		write(99,*)NAtoms
		write(99,*)
        do i=1,NAtoms
            write(99,'(A2,3F20.15)')MoleculeDetail.ElementSymbol(i),r(3*i-2:3*i)
        end do
    close(99)
    open(unit=99,file='MexInternalGeometry.out',status='replace')
        write(99,*)q
    close(99)
	open(unit=99,file='Mexg.out',status='replace')
	    write(99,*)(dH(:,InterestingState+1,InterestingState+1)-dH(:,InterestingState,InterestingState))/2d0
    close(99)
	open(unit=99,file='Mexh.out',status='replace')
	    write(99,*)dH(:,InterestingState+1,InterestingState)
    close(99)
    contains!Special routine for 2 state mex search
        subroutine f(Hd11,q,intdim)
            integer,intent(in)::intdim
            real*8,dimension(intdim),intent(in)::q
            real*8,intent(out)::Hd11
            real*8,dimension(NStates,NStates)::H
            H=Hd(q)
            Hd11=H(1,1)
        end subroutine f
        subroutine fd(dHd11,q,intdim)
            integer,intent(in)::intdim
            real*8,dimension(intdim),intent(in)::q
            real*8,dimension(intdim),intent(out)::dHd11
            real*8,dimension(intdim,NStates,NStates)::dH
            dH=dHd(q)
            dHd11=dH(:,1,1)
        end subroutine fd
        integer function fdd(ddHd11,q,intdim)
            integer,intent(in)::intdim
            real*8,dimension(intdim),intent(in)::q
            real*8,dimension(intdim,intdim),intent(out)::ddHd11
            real*8,dimension(intdim,intdim,NStates,NStates)::ddH
            ddH=ddHd(q)
            ddHd11=ddH(:,:,1,1)
            fdd=0!return 0
        end function fdd
        subroutine c(cx,q,M,intdim)
            integer,intent(in)::M,intdim
            real*8,dimension(intdim),intent(in)::q
            real*8,dimension(2),intent(out)::cx
            real*8,dimension(NStates,NStates)::H
            H=Hd(q)
            cx(1)=H(2,2)-H(1,1)
            cx(2)=H(2,1)
        end subroutine c
        subroutine cd(cdx,q,M,intdim)
            integer,intent(in)::M,intdim
            real*8,dimension(intdim),intent(in)::q
            real*8,dimension(intdim,2),intent(out)::cdx
            real*8,dimension(intdim,NStates,NStates)::dH
            dH=dHd(q)
            cdx(:,1)=dH(:,2,2)-dH(:,1,1)
            cdx(:,2)=dH(:,2,1)
        end subroutine cd
        integer function cdd(cddx,q,M,intdim)
            integer,intent(in)::M,intdim
            real*8,dimension(intdim),intent(in)::q
            real*8,dimension(intdim,intdim,2),intent(out)::cddx
            real*8,dimension(intdim,intdim,NStates,NStates)::ddH
            ddH=ddHd(q)
            cddx(:,:,1)=ddH(:,:,2,2)-ddH(:,:,1,1)
            cddx(:,:,2)=ddH(:,:,2,1)
            cdd=0!return 0
        end function cdd
end subroutine MexSearch

subroutine Evaluate()
end subroutine Evaluate

!---------- Auxiliary routine ----------
    !Reformat routine for calling nonlinear optimizers

    subroutine AdiabaticEnergyInterface(E,q,intdim)
        real*8,intent(out)::E
        integer,intent(in)::intdim
        real*8,dimension(intdim),intent(in)::q
        real*8,dimension(NStates)::energy
        energy=AdiabaticEnergy(q)
        E=energy(InterestingState)
    end subroutine AdiabaticEnergyInterface

    subroutine AdiabaticGradientInterface(dV,q,intdim)
        integer,intent(in)::intdim
        real*8,dimension(intdim),intent(out)::dV
        real*8,dimension(intdim),intent(in)::q
        real*8,dimension(intdim,NStates,NStates)::dH
        dH=AdiabaticdH(q)
        dV=dH(:,InterestingState,InterestingState)
    end subroutine AdiabaticGradientInterface

    integer function AdiabaticEnergy_GradientInterface(E,dV,q,intdim)
        real*8,intent(out)::E
        integer,intent(in)::intdim
        real*8,dimension(intdim),intent(out)::dV
        real*8,dimension(intdim),intent(in)::q
        real*8,dimension(NStates)::energy
        real*8,dimension(intdim,NStates,NStates)::dH
        call AdiabaticEnergy_dH(q,energy,dH)
        E=energy(InterestingState)
        dV=dH(:,InterestingState,InterestingState)
        AdiabaticEnergy_GradientInterface=0!return 0
    end function AdiabaticEnergy_GradientInterface
    
    integer function AdiabaticHessianInterface(Hessian,q,intdim)
        integer,intent(in)::intdim
        real*8,dimension(intdim,intdim),intent(out)::Hessian
        real*8,dimension(intdim),intent(in)::q
        real*8,dimension(intdim,intdim,NStates,NStates)::ddH
        ddH=AdiabaticddH(q)
        Hessian=ddH(:,:,InterestingState,InterestingState)
        AdiabaticHessianInterface=0!return 0
	end function AdiabaticHessianInterface
	
	subroutine AdiabaticGapInterface(gap,q,M,intdim)
		integer,intent(in)::M,intdim
		real*8,dimension(intdim),intent(in)::q
		real*8,dimension(1),intent(out)::gap
		real*8,dimension(NStates)::energy
		energy=AdiabaticEnergy(q)
		gap(1)=energy(InterestingState+1)-energy(InterestingState)
	end subroutine AdiabaticGapInterface

	subroutine AdiabaticGapGradientInterface(dgap,q,M,intdim)
		integer,intent(in)::M,intdim
		real*8,dimension(intdim),intent(in)::q
		real*8,dimension(intdim,1),intent(out)::dgap
		real*8,dimension(intdim,NStates,NStates)::dH
		dH=AdiabaticdH(q)
		dgap(:,1)=dH(:,InterestingState+1,InterestingState+1)-dH(:,InterestingState,InterestingState)
	end subroutine AdiabaticGapGradientInterface

	integer function AdiabaticGapHessianInterface(ddgap,q,M,intdim)
	    integer,intent(in)::M,intdim
		real*8,dimension(intdim),intent(in)::q
		real*8,dimension(intdim,intdim,1),intent(out)::ddgap
		real*8,dimension(intdim,intdim,NStates,NStates)::ddH
		ddH=AdiabaticddH(q)
		ddgap(:,:,1)=ddH(:,:,InterestingState+1,InterestingState+1)-ddH(:,:,InterestingState,InterestingState)
	    AdiabaticGapHessianInterface=0!return 0
    end function AdiabaticGapHessianInterface
!----------------- End -----------------

end module Analyzation