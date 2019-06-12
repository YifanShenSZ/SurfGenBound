!Provide analyzation of the fitted potential energy surface:
!    Minimum search and vibration analysis of specified adiabatic state
!    Mex search and gh orthogonalization along conical intersection seam between specified adiabatic states
!    Compute H & ▽H in diabatic & adiabatic & nondegenerate representation on specified geometries
module Analyzation
    use Basic
    implicit none

!Parameter
	real*8::Analyzation_NGrid=10,&!Generate how many grid points per direction
	        Analyzation_ghstep=0.01d0,&!Every how much bohr generate a grid point
			Analyzation_miu0=1d0!Initial strength of constraint violation penalty
	!Search control
		character*32::Analyzation_Searcher='BFGS'!Available: NewtonRaphson, BFGS, LBFGS, ConjugateGradient
	    logical::Analyzation_UseStrongWolfe=.true.!Whether use strong Wolfe condition instead of Wolfe condition
	    character*32::Analyzation_ConjugateGradientSolver='DY'!Available: DY (Dai-Yun), PR (Polak-Ribiere+)

!Analyzation module only variable
	!Modulewide accessed input variable
	    integer::Analyzation_state
	    logical::Analyzation_SearchDiabatic
	!Geometry information
		integer::Analyzation_NGeoms
		real*8,allocatable,dimension(:,:)::Analyzation_cartgeom,Analyzation_intgeom
		real*8,allocatable,dimension(:,:,:)::Analyzation_B
		real*8,allocatable,dimension(:)::Analyzation_g,Analyzation_h
    
contains
subroutine Analyze()!Top level standard interface for other modules to call
	character*32::JobType
	call ReadAnalyzeInput()
	select case(JobType)
	    case('evaluate')
	    	call evaluate()
        case('min')
            call MinimumSearch()
        case('mex')
            call MexSearch()
		case('OriginShift')
			call OriginShift(Analyzation_intgeom(:,1)-ReferencePoint.geom)
			JobType='HdNewOrigin.CheckPoint'
			call WriteHdExpansionCoefficients(JobType)
        case default!Throw a warning
            write(*,*)'Program abort: unsupported analyzation job type '//trim(adjustl(JobType))
            stop
	end select
	contains
	    subroutine ReadAnalyzeInput()!Read the input file for Analyze: analyzation.in
	    	logical::intgeom
	    	character*128::GeomFile,RefghFile,DispFile
	    	integer::i,j
	    	open(unit=99,file='analyzation.in',status='old')!Read main input, write some job comment
	    		read(99,*)
	    		read(99,*)
	    		read(99,*)
	    		read(99,*)JobType
	    			write(*,*)'Analyzation job type: '//JobType
	    		read(99,*)
	    		read(99,*)Analyzation_state
	    		read(99,*)
	    		read(99,*)GeomFile
	    		read(99,*)
				read(99,*)intgeom
				read(99,*)
	    		read(99,*)Analyzation_SearchDiabatic
	    		read(99,*)
	    		read(99,*)RefghFile
	    	close(99)
	    	open(unit=99,file=GeomFile,status='old')!Read geometries of interest and convert to internal coordinate
	    		Analyzation_NGeoms=0!Count number of geometries
	    		do
	    			read(99,*,iostat=i)
	    			if(i/=0) exit
	    			Analyzation_NGeoms=Analyzation_NGeoms+1
	    		end do
	    		if(intgeom) then
	    			Analyzation_NGeoms=Analyzation_NGeoms/InternalDimension
	    			rewind 99
	    			allocate(Analyzation_intgeom(InternalDimension,Analyzation_NGeoms))!Read geometries
	    			do i=1,Analyzation_NGeoms
	    				do j=1,InternalDimension
	    					read(99,*)Analyzation_intgeom(j,i)
	    				end do
	    			end do
	    		else
	    			Analyzation_NGeoms=Analyzation_NGeoms/MoleculeDetail.NAtoms
	    			rewind 99
	    			allocate(Analyzation_cartgeom(CartesianDimension,Analyzation_NGeoms))!Read geometries
	    			do i=1,Analyzation_NGeoms
	    				read(99,*)Analyzation_cartgeom(:,i)
	    			end do
	    			allocate(Analyzation_intgeom(InternalDimension,Analyzation_NGeoms))!Cart2int
	    			allocate(Analyzation_B(InternalDimension,CartesianDimension,Analyzation_NGeoms))
	    			do i=1,Analyzation_NGeoms
	    				call WilsonBMatrixAndInternalCoordinateq(Analyzation_B(:,:,i),Analyzation_intgeom(:,i),Analyzation_cartgeom(:,i),InternalDimension,CartesianDimension)
	    				Analyzation_intgeom(:,i)=Analyzation_intgeom(:,i)-ReferencePoint.geom!This program requires only internal coordinate difference
	    			end do
	    			open(unit=100,file='int'//trim(GeomFile)//'.out',status='replace')!Output an internal coordinate version for future use
	    				do i=1,Analyzation_NGeoms
	    					do j=1,InternalDimension
	    						write(100,*)Analyzation_intgeom(j,i)
	    					end do
	    				end do
	    			close(99)
	    		end if
	    	close(99)
	    	if(JobType=='mex') then!Look for reference g & h
	    		open(unit=99,file=RefghFile,status='old',iostat=i)
	    			if(i==0) then
	    				allocate(Analyzation_g(CartesianDimension))
	    				allocate(Analyzation_h(CartesianDimension))
	    				read(99,*)Analyzation_g
	    				read(99,*)Analyzation_h
	    			end if
	    		close(99)
	    	end if
	    end subroutine ReadAnalyzeInput
end subroutine Analyze

subroutine Evaluate()
	integer::i,j,k
	real*8,dimension(NState,Analyzation_NGeoms)::PES,ndSurface
	real*8,dimension(NState,NState,Analyzation_NGeoms)::HdSurface,dHdNormSurface,dHaNormSurface,dHndNormSurface
    real*8,dimension(NState,NState)::eigvec
	real*8,dimension(InternalDimension,NState,NState)::dH,dHa
	if(allocated(Analyzation_cartgeom)) then!Return Cartesian gradient
		do i=1,Analyzation_NGeoms!Compute
            HdSurface(:,:,i)=Hd(Analyzation_intgeom(:,i))
	    	dH=dHd(Analyzation_intgeom(:,i))
	    	forall(j=1:NState,k=1:NState,j>=k)
	    		dHdNormSurface(j,k,i)=norm2(matmul(transpose(Analyzation_B(:,:,i)),dH(:,j,k)))
	    	end forall
            eigvec=HdSurface(:,:,i)
            call My_dsyev('V',eigvec,PES(:,i),NState)
	    	dHa=sy3UnitaryTransformation(dH,eigvec,InternalDimension,NState)
	    	forall(j=1:NState,k=1:NState,j>=k)
	    	    dHaNormSurface(j,k,i)=norm2(matmul(transpose(Analyzation_B(:,:,i)),dHa(:,j,k)))
	    	end forall
	    	call NondegenerateRepresentation(dH,ndSurface(:,i),eigvec,InternalDimension,NState,DegenerateThreshold=AlmostDegenerate)
	    	forall(j=1:NState,k=1:NState,j>=k)
		        dHndNormSurface(j,k,i)=norm2(matmul(transpose(Analyzation_B(:,:,i)),dH(:,j,k)))
	    	end forall
		end do
	else!Return internal gradient
	    do i=1,Analyzation_NGeoms!Compute
            HdSurface(:,:,i)=Hd(Analyzation_intgeom(:,i))
	    	dH=dHd(Analyzation_intgeom(:,i))
	    	forall(j=1:NState,k=1:NState,j>=k)
	    		dHdNormSurface(j,k,i)=norm2(dH(:,j,k))
	    	end forall
            eigvec=HdSurface(:,:,i)
            call My_dsyev('V',eigvec,PES(:,i),NState)
	    	dHa=sy3UnitaryTransformation(dH,eigvec,InternalDimension,NState)
	    	forall(j=1:NState,k=1:NState,j>=k)
	    	    dHaNormSurface(j,k,i)=norm2(dHa(:,j,k))
	    	end forall
	    	call NondegenerateRepresentation(dH,ndSurface(:,i),eigvec,InternalDimension,NState,DegenerateThreshold=AlmostDegenerate)
	    	forall(j=1:NState,k=1:NState,j>=k)
		        dHndNormSurface(j,k,i)=norm2(dH(:,j,k))
	    	end forall
		end do
	end if
    !Output
    open(unit=99,file='PotentialEnergySurface.txt',status='replace')
        write(99,'(A9,A1)',advance='no')'Geometry#',char(9)
        do i=1,NState-1
            write(99,'(2x,A6,I2,A5,3x,A1)',advance='no')'Energy',i,'/cm-1',char(9)
        end do
        write(99,'(2x,A6,I2,A5,3x)')'Energy',NState,'/cm-1'
        do i=1,Analyzation_NGeoms
            write(99,'(I9,A1)',advance='no')i,char(9)
            do j=1,NState-1
                write(99,'(F18.8,A1)',advance='no')PES(j,i)/cm_1InAU,char(9)
            end do
            write(99,'(F18.8)')PES(j,i)/cm_1InAU
        end do
    close(99)
    open(unit=99,file='NondegenerateEigenValue.txt',status='replace')
        write(99,'(A9,A1)',advance='no')'Geometry#',char(9)
        do i=1,NState-1
            write(99,'(A10,I2,A5,1x,A1)',advance='no')'Eigenvalue',i,'/a.u.',char(9)
        end do
        write(99,'(A10,I2,A5,1x)')'Eigenvalue',i,'/a.u.'
        do i=1,Analyzation_NGeoms
            write(99,'(I9,A1)',advance='no')i,char(9)
            do j=1,NState-1
                write(99,'(F18.8,A1)',advance='no')ndSurface(j,i),char(9)
            end do
            write(99,'(F18.8)')ndSurface(j,i)
        end do
    close(99)
    open(unit=99,file='Hd.txt',status='replace')
        write(99,'(A10)',advance='no')'Geometry#'//char(9)
        do i=1,NState
            write(99,'(1x,A8,I2,A5,2x,A1)',advance='no')'Diagonal',i,'/cm-1',char(9)
		end do
		do i=1,NState
			do j=i+1,NState
				write(99,'(A8,I2,A1,I2,A5,A1)',advance='no')'Coupling',j,'&',i,'/cm-1',char(9)
			end do
        end do
        write(99,*)
        do i=1,Analyzation_NGeoms
            write(99,'(I9,A1)',advance='no')i,char(9)
            do j=1,NState
                write(99,'(F18.8,A1)',advance='no')HdSurface(j,j,i)/cm_1InAU,char(9)
			end do
			do k=1,NState
				do j=k+1,NState
					write(99,'(F18.8,A1)',advance='no')HdSurface(j,k,i)/cm_1InAU,char(9)
				end do
			end do
            write(99,*)
        end do
	close(99)
	open(unit=99,file='HdGradient.txt',status='replace')
        write(99,'(A10)',advance='no')'Geometry#'//char(9)
        do i=1,NState
            write(99,'(1x,A8,I2,A5,2x,A1)',advance='no')'Diagonal',i,'/a.u.',char(9)
		end do
		do i=1,NState
			do j=i+1,NState
				write(99,'(A8,I2,A1,I2,A5,A1)',advance='no')'Coupling',j,'&',i,'/a.u.',char(9)
			end do
        end do
        write(99,*)
        do i=1,Analyzation_NGeoms
            write(99,'(I9,A1)',advance='no')i,char(9)
            do j=1,NState
                write(99,'(F18.8,A1)',advance='no')dHdNormSurface(j,j,i),char(9)
			end do
			do k=1,NState
				do j=k+1,NState
					write(99,'(F18.8,A1)',advance='no')dHdNormSurface(j,k,i),char(9)
				end do
			end do
            write(99,*)
        end do
	close(99)
	open(unit=99,file='HaGradient.txt',status='replace')
        write(99,'(A10)',advance='no')'Geometry#'//char(9)
        do i=1,NState
            write(99,'(A11,I2,A5,A1)',advance='no')'Energy grad',i,'/a.u.',char(9)
		end do
		do i=1,NState
			do j=i+1,NState
				write(99,'(2x,A3,I2,A1,I2,A5,3x,A1)',advance='no')'ISC',j,'&',i,'/a.u.',char(9)
			end do
        end do
		do i=1,NState
			do j=i+1,NState
				write(99,'(2x,A3,I2,A1,I2,A5,3x,A1)',advance='no')'NAC',j,'&',i,'/a.u.',char(9)
			end do
		end do
        write(99,*)
        do i=1,Analyzation_NGeoms
            write(99,'(I9,A1)',advance='no')i,char(9)
            do j=1,NState
                write(99,'(F18.8,A1)',advance='no')dHaNormSurface(j,j,i),char(9)
			end do
			do k=1,NState
				do j=k+1,NState
					write(99,'(F18.8,A1)',advance='no')dHaNormSurface(j,k,i),char(9)
				end do
			end do
			do k=1,NState
				do j=k+1,NState
					write(99,'(F18.8,A1)',advance='no')dHaNormSurface(j,k,i)/dABS(PES(k,i)-PES(j,i)),char(9)
				end do
			end do
            write(99,*)
        end do
	close(99)
	open(unit=99,file='HndGradient.txt',status='replace')
        write(99,'(A10)',advance='no')'Geometry#'//char(9)
        do i=1,NState
            write(99,'(1x,A8,I2,A5,2x,A1)',advance='no')'Diagonal',i,'/a.u.',char(9)
		end do
		do i=1,NState
			do j=i+1,NState
				write(99,'(A8,I2,A1,I2,A5,A1)',advance='no')'Coupling',j,'&',i,'/a.u.',char(9)
			end do
        end do
        write(99,*)
        do i=1,Analyzation_NGeoms
            write(99,'(I9,A1)',advance='no')i,char(9)
            do j=1,NState
                write(99,'(F18.8,A1)',advance='no')dHndNormSurface(j,j,i),char(9)
			end do
			do k=1,NState
				do j=k+1,NState
					write(99,'(F18.8,A1)',advance='no')dHndNormSurface(j,k,i)/dABS(ndSurface(k,i)-ndSurface(j,i)),char(9)
				end do
			end do
            write(99,*)
        end do
	close(99)
end subroutine Evaluate

subroutine MinimumSearch()
    integer::i,j
    real*8,dimension(InternalDimension)::q,freq
    real*8,dimension(CartesianDimension)::r,rtemp
    real*8,dimension(InternalDimension,CartesianDimension)::B
	real*8,dimension(InternalDimension,InternalDimension)::Hessian,mode
	write(*,'(1x,A46,1x,I2)')'Search for minimum on potential energy surface',Analyzation_state
	q=Analyzation_intgeom(:,1)
	if(Analyzation_SearchDiabatic) then
	    select case(Analyzation_Searcher)
        	case('NewtonRaphson')
        		call NewtonRaphson(DiabaticEnergyInterface,DiabaticGradientInterface,q,InternalDimension,&
        			fdd=DiabaticHessianInterface,Strong=Analyzation_UseStrongWolfe)
        	case('BFGS')
        		call BFGS(DiabaticEnergyInterface,DiabaticGradientInterface,q,InternalDimension,&
        			fdd=DiabaticHessianInterface,Strong=Analyzation_UseStrongWolfe)
        	case('LBFGS')
        		call LBFGS(DiabaticEnergyInterface,DiabaticGradientInterface,q,InternalDimension,&
        			Strong=Analyzation_UseStrongWolfe)
        	case('ConjugateGradient')
        		call ConjugateGradient(DiabaticEnergyInterface,DiabaticGradientInterface,q,InternalDimension,&
        			Strong=Analyzation_UseStrongWolfe,Method=Analyzation_ConjugateGradientSolver)
        	case default!Throw a warning
        		write(*,*)'Program abort: unsupported searcher '//trim(adjustl(Analyzation_Searcher))
        		stop
        end select
	else
	    select case(Analyzation_Searcher)
	        case('NewtonRaphson')
	        	call NewtonRaphson(AdiabaticEnergyInterface,AdiabaticGradientInterface,q,InternalDimension,&
	        		fdd=AdiabaticHessianInterface,f_fd=AdiabaticEnergy_GradientInterface,&
	        		Strong=Analyzation_UseStrongWolfe)
            case('BFGS')
                call BFGS(AdiabaticEnergyInterface,AdiabaticGradientInterface,q,InternalDimension,&
	        		fdd=AdiabaticHessianInterface,f_fd=AdiabaticEnergy_GradientInterface,&
	        		Strong=Analyzation_UseStrongWolfe)
	        case('LBFGS')
	        	call LBFGS(AdiabaticEnergyInterface,AdiabaticGradientInterface,q,InternalDimension,&
	        		f_fd=AdiabaticEnergy_GradientInterface,Strong=Analyzation_UseStrongWolfe)
	        case('ConjugateGradient')
	        	call ConjugateGradient(AdiabaticEnergyInterface,AdiabaticGradientInterface,q,InternalDimension,&
	        		f_fd=AdiabaticEnergy_GradientInterface,Strong=Analyzation_UseStrongWolfe,&
	        		Method=Analyzation_ConjugateGradientSolver)
	        case default!Throw a warning
	        	write(*,*)'Program abort: unsupported searcher '//trim(adjustl(Analyzation_Searcher))
	        	stop
		end select
	end if
	open(unit=99,file='MinimumInternalGeometry.out',status='replace')
        write(99,*)q
	close(99)
	i=AdiabaticHessianInterface(Hessian,q,InternalDimension)
	q=q+ReferencePoint.geom
	if(allocated(Analyzation_cartgeom)) then
		rtemp=Analyzation_cartgeom(:,1)
	else
		rtemp=reshape(MoleculeDetail.RefConfig,[CartesianDimension])
	end if
	call StandardizeGeometry(rtemp,MoleculeDetail.mass,MoleculeDetail.NAtoms,1)
	r=CartesianCoordinater(q,CartesianDimension,InternalDimension,&
		mass=MoleculeDetail.mass,r0=rtemp)
	open(unit=99,file='MinimumCartesianGeometry.xyz',status='replace')
		write(99,*)MoleculeDetail.NAtoms
		write(99,*)
        do i=1,MoleculeDetail.NAtoms
            write(99,'(A2,3F20.15)')MoleculeDetail.ElementSymbol(i),r(3*i-2:3*i)/AInAU
        end do
    close(99)
    call WilsonBMatrixAndInternalCoordinateq(B,q,r,InternalDimension,CartesianDimension)
    call WilsonGFMethod(freq,mode,Hessian,InternalDimension,B,MoleculeDetail.mass,MoleculeDetail.NAtoms)
    open(unit=99,file='VibrationalFrequency.txt',status='replace')
        write(99,'(A4,A1,A14)')'Mode',char(9),'Frequency/cm-1'
        do i=1,InternalDimension
            write(99,'(I4,A1,F14.8)')i,char(9),freq(i)/cm_1InAu
        end do
    close(99)
	open(unit=99,file='NormalMode.txt',status='replace')
	    write(99,'(A6,A1)',advance='no')'q\Mode',char(9)
		do i=1,InternalDimension-1
			write(99,'(I6,A1)',advance='no')i,char(9)
		end do
		write(99,'(I6)')InternalDimension
		do i=1,InternalDimension
			write(99,'(I8,A1)',advance='no')i,char(9)
			do j=1,InternalDimension-1
				write(99,'(F18.15,A1)',advance='no')mode(j,i),char(9)
			end do
			write(99,'(F18.15)')mode(InternalDimension,i)
		end do
    close(99)
end subroutine MinimumSearch

subroutine MexSearch()
	integer::i,j
    real*8,dimension(InternalDimension)::q
	real*8,dimension(CartesianDimension)::r,rtemp,g,h
	real*8,dimension(InternalDimension,NState,NState)::intdH
	real*8,dimension(CartesianDimension,NState,NState)::cartdH
	if(allocated(Analyzation_g).and.allocated(Analyzation_h)) then!gh path for Columbus
		r=Analyzation_cartgeom(:,1)
		g=Analyzation_g/norm2(Analyzation_g)
		h=Analyzation_h/norm2(Analyzation_h)
		open(unit=99,file='gPath.geom',status='replace')
		    do i=-Analyzation_NGrid,-1
	        	rtemp=r+dble(i)*Analyzation_ghstep*g
		    	do j=1,MoleculeDetail.NAtoms
		    		write(99,'(A2,I8,3F14.8,F14.8)')MoleculeDetail.ElementSymbol(j),Symbol2Number(MoleculeDetail.ElementSymbol(j)),rtemp(3*j-2:3*j),MoleculeDetail.mass(j)/AMUInAU
		    	end do
			end do
			do i=1,Analyzation_NGrid
	        	rtemp=r+dble(i)*Analyzation_ghstep*g
		    	do j=1,MoleculeDetail.NAtoms
		    		write(99,'(A2,I8,3F14.8,F14.8)')MoleculeDetail.ElementSymbol(j),Symbol2Number(MoleculeDetail.ElementSymbol(j)),rtemp(3*j-2:3*j),MoleculeDetail.mass(j)/AMUInAU
		    	end do
		    end do
		close(99)
		open(unit=99,file='hPath.geom',status='replace')
		    do i=-Analyzation_NGrid,-1
	        	rtemp=r+dble(i)*Analyzation_ghstep*h
		    	do j=1,MoleculeDetail.NAtoms
		    		write(99,'(A2,I8,3F14.8,F14.8)')MoleculeDetail.ElementSymbol(j),Symbol2Number(MoleculeDetail.ElementSymbol(j)),rtemp(3*j-2:3*j),MoleculeDetail.mass(j)/AMUInAU
		    	end do
			end do
			do i=1,Analyzation_NGrid
	        	rtemp=r+dble(i)*Analyzation_ghstep*h
		    	do j=1,MoleculeDetail.NAtoms
		    		write(99,'(A2,I8,3F14.8,F14.8)')MoleculeDetail.ElementSymbol(j),Symbol2Number(MoleculeDetail.ElementSymbol(j)),rtemp(3*j-2:3*j),MoleculeDetail.mass(j)/AMUInAU
		    	end do
		    end do
		close(99)
	end if
	write(*,'(1x,A48,1x,I2,1x,A3,1x,I2)')'Search for mex between potential energy surfaces',Analyzation_state,'and',Analyzation_state+1
    q=Analyzation_intgeom(:,1)
	if(NState==2.and.Analyzation_SearchDiabatic) then!2 state case we can simply search for minimum of Hd diagonal subject to zero off-diagonal and degenerate diagonals
		call AugmentedLagrangian(f,fd,c,cd,q,InternalDimension,2,fdd=fdd,cdd=cdd,Precision=1d-8,&!This is Columbus7 energy precision
		    miu0=Analyzation_miu0,UnconstrainedSolver=Analyzation_Searcher,Method=Analyzation_ConjugateGradientSolver)
    else!In general case we have to search for minimum on potential energy surface of interest subject to degeneracy constaint
	    call AugmentedLagrangian(AdiabaticEnergyInterface,AdiabaticGradientInterface,AdiabaticGapInterface,AdiabaticGapGradientInterface,&
	        q,InternalDimension,1,Precision=1d-8,&!This is Columbus7 energy precision
			fdd=AdiabaticHessianInterface,cdd=AdiabaticGapHessianInterface,f_fd=AdiabaticEnergy_GradientInterface,&
			miu0=Analyzation_miu0,UnconstrainedSolver=Analyzation_Searcher,Method=Analyzation_ConjugateGradientSolver)
	end if
	open(unit=99,file='MexInternalGeometry.out',status='replace')
        write(99,*)q
    close(99)
	intdH=AdiabaticdH(q)
	q=q+ReferencePoint.geom
	if(allocated(Analyzation_cartgeom)) then
		rtemp=Analyzation_cartgeom(:,1)
	else
		rtemp=reshape(MoleculeDetail.RefConfig,[CartesianDimension])
	end if
	call StandardizeGeometry(rtemp,MoleculeDetail.mass,MoleculeDetail.NAtoms,1)
	call Internal2Cartesian(q,InternalDimension,r,CartesianDimension,NState,&
	    intnadgrad=intdH,cartnadgrad=cartdH,&
		mass=MoleculeDetail.mass,r0=rtemp)
	if(allocated(Analyzation_g).and.allocated(Analyzation_h)) then
		call ghOrthogonalization(cartdH(:,Analyzation_state,Analyzation_state),cartdH(:,Analyzation_state+1,Analyzation_state+1),cartdH(:,Analyzation_state+1,Analyzation_state),CartesianDimension,&
			gref=Analyzation_g,href=Analyzation_h)
	else
		call ghOrthogonalization(cartdH(:,Analyzation_state,Analyzation_state),cartdH(:,Analyzation_state+1,Analyzation_state+1),cartdH(:,Analyzation_state+1,Analyzation_state),CartesianDimension)
	end if
	g=(cartdH(:,Analyzation_state+1,Analyzation_state+1)-cartdH(:,Analyzation_state,Analyzation_state))/2d0
	h=cartdH(:,Analyzation_state+1,Analyzation_state)
	open(unit=99,file='MexCartesianGeometry.xyz',status='replace')
		write(99,*)MoleculeDetail.NAtoms
		write(99,*)
        do i=1,MoleculeDetail.NAtoms
            write(99,'(A2,3F20.15)')MoleculeDetail.ElementSymbol(i),r(3*i-2:3*i)/AInAU
        end do
	close(99)
	open(unit=99,file='Mexg.out',status='replace')
	    write(99,*)g
    close(99)
	open(unit=99,file='Mexh.out',status='replace')
	    write(99,*)h
	close(99)
	g=g/norm2(g)
	h=h/norm2(h)
	open(unit=99,file='gPathToEvaluate.in',status='replace')
	    do i=-Analyzation_NGrid,Analyzation_NGrid
	    	rtemp=r+dble(i)*Analyzation_ghstep*g
			write(99,*)rtemp
		end do
	close(99)
	open(unit=99,file='hPathToEvaluate.in',status='replace')
	    do i=-Analyzation_NGrid,Analyzation_NGrid
	    	rtemp=r+dble(i)*Analyzation_ghstep*h
			write(99,*)rtemp
		end do
	close(99)
	open(unit=99,file='DoubleConeToEvaluate.in',status='replace')
		do i=-Analyzation_NGrid,Analyzation_NGrid
			do j=-Analyzation_NGrid,Analyzation_NGrid
				rtemp=r+dble(i)*Analyzation_ghstep*g+dble(j)*Analyzation_ghstep*h
			    write(99,*)rtemp
			end do
		end do
	close(99)
    contains!Special routine for 2 state mex search
        subroutine f(Hd11,q,intdim)
            integer,intent(in)::intdim
            real*8,dimension(intdim),intent(in)::q
            real*8,intent(out)::Hd11
            real*8,dimension(NState,NState)::H
            H=Hd(q)
            Hd11=H(1,1)
        end subroutine f
        subroutine fd(dHd11,q,intdim)
            integer,intent(in)::intdim
            real*8,dimension(intdim),intent(in)::q
            real*8,dimension(intdim),intent(out)::dHd11
            real*8,dimension(intdim,NState,NState)::dH
            dH=dHd(q)
            dHd11=dH(:,1,1)
        end subroutine fd
        integer function fdd(ddHd11,q,intdim)
            integer,intent(in)::intdim
            real*8,dimension(intdim),intent(in)::q
            real*8,dimension(intdim,intdim),intent(out)::ddHd11
            real*8,dimension(intdim,intdim,NState,NState)::ddH
            ddH=ddHd(q)
            ddHd11=ddH(:,:,1,1)
            fdd=0!return 0
        end function fdd
        subroutine c(cx,q,M,intdim)
            integer,intent(in)::M,intdim
            real*8,dimension(intdim),intent(in)::q
            real*8,dimension(2),intent(out)::cx
            real*8,dimension(NState,NState)::H
            H=Hd(q)
            cx(1)=H(2,2)-H(1,1)
            cx(2)=H(2,1)
        end subroutine c
        subroutine cd(cdx,q,M,intdim)
            integer,intent(in)::M,intdim
            real*8,dimension(intdim),intent(in)::q
            real*8,dimension(intdim,2),intent(out)::cdx
            real*8,dimension(intdim,NState,NState)::dH
            dH=dHd(q)
            cdx(:,1)=dH(:,2,2)-dH(:,1,1)
            cdx(:,2)=dH(:,2,1)
        end subroutine cd
        integer function cdd(cddx,q,M,intdim)
            integer,intent(in)::M,intdim
            real*8,dimension(intdim),intent(in)::q
            real*8,dimension(intdim,intdim,2),intent(out)::cddx
            real*8,dimension(intdim,intdim,NState,NState)::ddH
            ddH=ddHd(q)
            cddx(:,:,1)=ddH(:,:,2,2)-ddH(:,:,1,1)
            cddx(:,:,2)=ddH(:,:,2,1)
            cdd=0!return 0
        end function cdd
end subroutine MexSearch

!---------- Auxiliary routine ----------
    !Reformat routine for calling nonlinear optimizers

    subroutine AdiabaticEnergyInterface(E,q,intdim)
        real*8,intent(out)::E
        integer,intent(in)::intdim
        real*8,dimension(intdim),intent(in)::q
        real*8,dimension(NState)::energy
        energy=AdiabaticEnergy(q)
        E=energy(Analyzation_state)
    end subroutine AdiabaticEnergyInterface

    subroutine AdiabaticGradientInterface(dV,q,intdim)
        integer,intent(in)::intdim
        real*8,dimension(intdim),intent(out)::dV
        real*8,dimension(intdim),intent(in)::q
        real*8,dimension(intdim,NState,NState)::dH
        dH=AdiabaticdH(q)
        dV=dH(:,Analyzation_state,Analyzation_state)
    end subroutine AdiabaticGradientInterface

    integer function AdiabaticEnergy_GradientInterface(E,dV,q,intdim)
        real*8,intent(out)::E
        integer,intent(in)::intdim
        real*8,dimension(intdim),intent(out)::dV
        real*8,dimension(intdim),intent(in)::q
        real*8,dimension(NState)::energy
        real*8,dimension(intdim,NState,NState)::dH
        call AdiabaticEnergy_dH(q,energy,dH)
        E=energy(Analyzation_state)
        dV=dH(:,Analyzation_state,Analyzation_state)
        AdiabaticEnergy_GradientInterface=0!return 0
    end function AdiabaticEnergy_GradientInterface
    
    integer function AdiabaticHessianInterface(Hessian,q,intdim)
        integer,intent(in)::intdim
        real*8,dimension(intdim,intdim),intent(out)::Hessian
        real*8,dimension(intdim),intent(in)::q
        real*8,dimension(intdim,intdim,NState,NState)::ddH
        ddH=AdiabaticddH(q)
        Hessian=ddH(:,:,Analyzation_state,Analyzation_state)
        AdiabaticHessianInterface=0!return 0
	end function AdiabaticHessianInterface
	
	subroutine AdiabaticGapInterface(gap,q,M,intdim)
		integer,intent(in)::M,intdim
		real*8,dimension(intdim),intent(in)::q
		real*8,dimension(1),intent(out)::gap
		real*8,dimension(NState)::energy
		energy=AdiabaticEnergy(q)
		gap(1)=energy(Analyzation_state+1)-energy(Analyzation_state)
	end subroutine AdiabaticGapInterface

	subroutine AdiabaticGapGradientInterface(dgap,q,M,intdim)
		integer,intent(in)::M,intdim
		real*8,dimension(intdim),intent(in)::q
		real*8,dimension(intdim,1),intent(out)::dgap
		real*8,dimension(intdim,NState,NState)::dH
		dH=AdiabaticdH(q)
		dgap(:,1)=dH(:,Analyzation_state+1,Analyzation_state+1)-dH(:,Analyzation_state,Analyzation_state)
	end subroutine AdiabaticGapGradientInterface

	integer function AdiabaticGapHessianInterface(ddgap,q,M,intdim)
	    integer,intent(in)::M,intdim
		real*8,dimension(intdim),intent(in)::q
		real*8,dimension(intdim,intdim,1),intent(out)::ddgap
		real*8,dimension(intdim,intdim,NState,NState)::ddH
		ddH=AdiabaticddH(q)
		ddgap(:,:,1)=ddH(:,:,Analyzation_state+1,Analyzation_state+1)-ddH(:,:,Analyzation_state,Analyzation_state)
	    AdiabaticGapHessianInterface=0!return 0
	end function AdiabaticGapHessianInterface
	
	subroutine DiabaticEnergyInterface(E,q,intdim)
		real*8,intent(out)::E
		integer,intent(in)::intdim
		real*8,dimension(intdim),intent(in)::q
		real*8,dimension(NState,NState)::H
		H=Hd(q)
		E=H(Analyzation_state,Analyzation_state)
	end subroutine DiabaticEnergyInterface

	subroutine DiabaticGradientInterface(dV,q,intdim)
        integer,intent(in)::intdim
        real*8,dimension(intdim),intent(out)::dV
        real*8,dimension(intdim),intent(in)::q
        real*8,dimension(intdim,NState,NState)::dH
        dH=dHd(q)
        dV=dH(:,Analyzation_state,Analyzation_state)
	end subroutine DiabaticGradientInterface
	
	integer function DiabaticHessianInterface(Hessian,q,intdim)
        integer,intent(in)::intdim
        real*8,dimension(intdim,intdim),intent(out)::Hessian
        real*8,dimension(intdim),intent(in)::q
        real*8,dimension(intdim,intdim,NState,NState)::ddH
        ddH=ddHd(q)
        Hessian=ddH(:,:,Analyzation_state,Analyzation_state)
        DiabaticHessianInterface=0!return 0
	end function DiabaticHessianInterface
!----------------- End -----------------

end module Analyzation