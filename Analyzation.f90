!Provide analyzation of the fitted potential energy surface:
!    Compute H & ▽H in diabatic & adiabatic & nondegenerate representation on specified geometries
!    Minimum search and vibration analysis of specified adiabatic state
!    Mex search and gh orthogonalization along conical intersection seam between specified adiabatic states
module Analyzation
    use Basic
    use DiabaticHamiltonian
    implicit none

!Analyzation module only variable
	!Analyzation input
        character*32::AnalyzationJobType
		integer::InterestingState
	!Geometry input
		integer::NGeoms
        real*8,allocatable,dimension(:,:)::InterestingGeom
    
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

subroutine ReadAnalyzeInput()!Read the input file for Analyzation: AnalyzeInput
	character*128::InterestingGeomFile
	integer::i
	real*8,dimension(3)::buffer
    open(unit=99,file='analyzation.in',status='old')
        read(99,*)
        read(99,*)
        read(99,*)
        read(99,*)AnalyzationJobType
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
    q=InternalCoordinateq(InterestingGeom(:,1),InternalDImension,CartesianDimension)
    call BFGS(AdiabaticEnergyInterface,AdiabaticGradientInterface,q,InternalDImension,&
        fdd=AdiabaticHessianInterface,f_fd=AdiabaticEnergy_GradientInterface,Strong=.true.)
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
    open(unit=99,file='MinimumNormalCoordinate.out',status='replace')
        write(99,*)Hessian
    close(99)
end subroutine MinimumSearch

subroutine MexSearch()

end subroutine MexSearch

subroutine Evaluate()
end subroutine Evaluate

!---------- Auxiliary routine ----------
    !Reformat routine for NonlinearOptimization

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
!----------------- End -----------------

end module Analyzation