!SurfGenBound: surface generation package for bounded system
!Construct diabatic hamiltonian (Hd) by least square fitting H & ▽H
!IO unit: atomic unit unless specified in file
!Computation unit: atomic unit
!              |- NadVibSInterface <-|<-     Analyzation     <-|
!Level: Main <-|-                    |<- ElectronicStructure <-|- Basic <- DiabaticHamiltonian
!              |-              HdLeastSquareFit              <-|
program main
    use Basic
    use HdLeastSquareFit; use ElectronicStructure; use Analyzation
    use NadVibSInterface
    implicit none
    !Main only accessed input variable
        character*32::JobType
        logical::SameTrainingSet
        integer::IndexReference
        character*128::ArtifactGeometryDataFile,ArtifactEnergyDataFile
    !Job control
        !Whether perform the fitting procedure gradually:
        !    1. Start from a several data points nearest to the reference point.
        !       The number points is either automatically determined or manually entred,
        !       and the program will guaranteen there are more equations than variables.
        !    2. Then the nearest of the remaining points will be added in. Repeat until all points are fitted.
        logical::GradualFit=.false.,AutoGradualFit=.true.
        integer::ManualNPoints
        !If gradual fit enabled, then if the data points are unaligned, align all geometries such that:
        !    the centre of mass is at origin
        !    the rotational principle axes are along xyz axes, with the smallest corresponding to x axis, 2nd to y, 3rd to z
        !        Since the positive direction cannot be determined yet, the final geometry will be the legal one
        !        with smallest difference to the reference geometry (difference = 2 norm square of Cartesian coordinate difference)
        !    the difference value is also the square of the Cartesian distance to the reference geometry
        logical::Unaligned=.true.
        !If gradual fit enabled, then if the data points are unsorted, sort according to the Cartesian distance to the reference geometry in ascending order
        !    Centre of mass position and orientation should not matter, so the points have to be aligned
        logical::Unsorted=.true.
    !Gradual fit variables
        integer::NPointsInput,NDegeneratePointsInput,NArtifactPointsInput!Store the input value (before change it)
		real*8::HdLSF_RegularizationOld!Store the original parameter value (before change it)
!---------- Initialize ----------
    write(*,'(A61)')'SurfGenBound: a surface generation package for bounded system'
    write(*,'(A15)')'Yifan Shen 2019'
    write(*,*)
    call ShowTime()
    write(*,*)
    write(*,*)'Electronic structure software in use = '//ElectronicStructureSoftware
    call ReadInput(); call Initialize()
!------------- End --------------

!----------- Run job ------------
    select case(JobType)
	case('FitNewDiabaticHamiltonian')
        if(GradualFit) then
            write(*,'(1x,A39)')'Perform the fitting procedure gradually'
            !Store the original values then change them
                NPointsInput=NPoints
                NDegeneratePointsInput=NDegeneratePoints
                    NDegeneratePoints=0d0!Almost degenerate points are omitted during the gradual fit procedure
                NArtifactPointsInput=NArtifactPoints
                    NArtifactPoints=0!Unreliable points are omitted during the gradual fit procedure
                HdLSF_RegularizationOld=HdLSF_Regularization
                    HdLSF_Regularization=0d0!Regularization is disabled during the gradual fit procedure
            !At least this number of data points are required to have more equations than variables
            NPoints=ceiling(dble(NHdExpansionCoefficients)/dble(DataPerPoint))
            if(AutoGradualFit) then
                write(*,'(1x,A54)')'Generating initial guess with fewest nearest points...'
            else
                if(ManualNPoints<NPoints) then
                    write(*,'(1x,A67)')'Warning: Too few data points, equations must be more than variables'
                else
                    NPoints=ManualNPoints
                end if
                write(*,'(1x,A72)')'Generating initial guess with user specified number of nearest points...'
            end if
            do while(NPoints<NPointsInput)
                write(*,*)
                write(*,'(1x,A33,1x,I6)')'Number of data points now in use:',NPoints
                call InitializeHdLeastSquareFit()
                call FitHd()
                write(*,'(1x,A25)')'Save guadual fit progress'
                open(unit=99,file='GradualFit.CheckPoint',status='replace')
                    write(99,'(A48)')'Number of data points to be least square fitted:'
                    write(99,*)NPoints
                close(99)
                NPoints=NPoints+1
            end do
            write(*,*)
            write(*,'(1x,A52)')'All data points added in, fitting Hd in usual way...'
                NDegeneratePoints=NDegeneratePointsInput
                NArtifactPoints=NArtifactPointsInput
                HdLSF_Regularization=HdLSF_RegularizationOld
            call InitializeHdLeastSquareFit()
            call FitHd()
        else
            write(*,*)'Fitting Hd...'
            call FitHd()
        end if
    case('ContinueFitting')
        if(GradualFit) then
            !Store the original values then change them
                NPointsInput=NPoints
                NDegeneratePointsInput=NDegeneratePoints
                    NDegeneratePoints=0d0!Almost degenerate points are omitted during the gradual fit procedure
                NArtifactPointsInput=NArtifactPoints
                    NArtifactPoints=0!Unreliable points are omitted during the gradual fit procedure
                HdLSF_RegularizationOld=HdLSF_Regularization
                    HdLSF_Regularization=0d0!Regularization is disabled during the gradual fit procedure
            write(*,'(1x,A33)')'Continue guadual fit procedure...'
            open(unit=99,file='GradualFit.CheckPoint',status='old')
                read(99,*)
                read(99,*)NPoints
            close(99)
            do while(NPoints<NPointsInput)
                write(*,*)
                write(*,'(1x,A33,1x,I6)')'Number of data points now in use:',NPoints
                call InitializeHdLeastSquareFit()
                call FitHd()
                write(*,'(1x,A25)')'Save guadual fit progress'
                open(unit=99,file='GradualFit.CheckPoint',status='replace')
                    write(99,'(A48)')'Number of data points to be least square fitted:'
                    write(99,*)NPoints
                close(99)
                NPoints=NPoints+1
            end do
            write(*,*)
            write(*,'(1x,A52)')'All data points added in, fitting Hd in usual way...'
                NDegeneratePoints=NDegeneratePointsInput
                NArtifactPoints=NArtifactPointsInput
                HdLSF_Regularization=HdLSF_RegularizationOld
            call InitializeHdLeastSquareFit()
            call FitHd()
        else
            write(*,*)'Fitting Hd...'
            call FitHd()
        end if
    case('Analyze'); call Analyze()
    case('NadVibS'); call GenerateNadVibSInput()
    case default; write(*,*)'Program abort: unsupported job type '//JobType; stop
    end select
!------------- End --------------

!---------- Clean up ------------
    call ShowTime(); write(*,'(1x,A16)')'Mission complete'
!------------- End --------------

contains
subroutine ReadInput()!Read main input files: SurfGenBound.in, eg.xyz, advance.in (optional)
    character*128::MoleculeDetailFile
    logical::advance
    integer::i
    open(unit=99,file='SurfGenBound.in',status='old')!Read main input, write some job comment
        read(99,*); read(99,*); read(99,*); read(99,*)JobType; write(*,*)'Job type: '//JobType
        read(99,*); read(99,*)MoleculeDetailFile
        read(99,*); read(99,*)NState
        read(99,*); read(99,*)SameTrainingSet
        read(99,*); read(99,*)NPoints
        read(99,*); read(99,*)IndexReference; if(IndexReference>NPoints) stop 'Program abort: the index of reference point must <= number of data points'
        read(99,*); read(99,*)ArtifactGeometryDataFile
        read(99,*); read(99,*)ArtifactEnergyDataFile
        read(99,*); read(99,*)NArtifactPoints
        read(99,*); read(99,*)advance
    close(99)
    open(unit=99,file=MoleculeDetailFile,status='old')!Read molecule detail
        read(99,*)MoleculeDetail.NAtoms
            allocate(MoleculeDetail.ElementSymbol(MoleculeDetail.NAtoms))
            allocate(MoleculeDetail.RefConfig(3,MoleculeDetail.NAtoms))
            allocate(MoleculeDetail.mass(MoleculeDetail.NAtoms))
        read(99,*)
        do i=1,MoleculeDetail.NAtoms
            read(99,*)MoleculeDetail.ElementSymbol(i),MoleculeDetail.RefConfig(:,i)
            MoleculeDetail.ElementSymbol(i)=trim(adjustl(MoleculeDetail.ElementSymbol(i)))
        end do
        MoleculeDetail.RefConfig=MoleculeDetail.RefConfig*AInAU!Convert to atomic unit
        read(99,*)
        do i=1,MoleculeDetail.NAtoms
            read(99,*)MoleculeDetail.mass(i)
        end do
        MoleculeDetail.mass=MoleculeDetail.mass*AMUInAU!Convert to atomic unit
    close(99)
    if(advance) then!If requested, read advanced input
        write(*,'(1x,A66)')'Advanced input requested, parameters are set to user specification'
        open(unit=99,file='advance.in',status='old')
            namelist /AdvancedInput/ &
                !Basic
                    HighEnergy,AlmostDegenerate,&
                !HdLeastSquareFit
                    HdLSF_Regularization,HdLSF_Solver,&
                    HdLSF_MaxHopperIteration,HdLSF_MaxLocalMinimizerIteration,HdLSF_Max2StepIteration,&
                    HdLSF_pseudolinearFollowStep,HdLSF_pseudolinearMaxMonotonicalIncrease,&
                    HdLSF_LineSearcher,HdLSF_UseStrongWolfe,HdLSF_LBFGSMemory,HdLSF_ConjugateGradientSolver,&
                !Analyzation
                    Analyzation_NGrid,Analyzation_ghstep,Analyzation_miu0,&
                    Analyzation_Searcher,Analyzation_UseStrongWolfe,Analyzation_ConjugateGradientSolver,&
                !NadVibSInterface
                    NVS_contour,&
                !Main
                    GradualFit,AutoGradualFit,ManualNPoints,Unaligned,Unsorted
            read(99,nml=AdvancedInput)
        close(99)
    end if
end subroutine ReadInput

subroutine Initialize()!Program initializer
	character*128::CharTemp128; logical::flag
    integer::istate,jstate,i,j; integer,dimension(1)::indice; real*8::dbletemp
    real*8,allocatable,dimension(:)::OldRefGeom,OldRefEnergy
    !General initialize
        call BetterRandomSeed()
        call InitializeBasic()
        CartesianDimension=3*MoleculeDetail.NAtoms
        InternalDimension=DefineInternalCoordinate(ElectronicStructureSoftware)
    select case(JobType)!Job specific initialize
    case('FitNewDiabaticHamiltonian')!To fit Hd from scratch, read training set then rearrange it
        call Initialize_NewTrainingSet()
		call InitializeDiabaticHamiltonian(NState,InternalDimension,NewHd=.true.)
		!Provide an initial guess of Hd
            call CheckDegeneracy(flag,AlmostDegenerate,ReferencePoint.energy,NState)
            i=WhichExpansionBasis(0,indice(1:0))
            if(i>0) then
                if(flag) then
                    forall(istate=1:NState,jstate=1:NState,istate>=jstate)
                        Hd_HdEC(istate,jstate).Array(i)=ReferencePoint.H(istate,jstate)
                    end forall
                else
                    forall(istate=2:NState)
                        Hd_HdEC(istate,istate).Array(i)=ReferencePoint.energy(istate)-ReferencePoint.energy(1)
                    end forall
                end if
            end if
            do jstate=1,NState
                do istate=jstate,NState
                    do j=1,InternalDimension
                        indice(1)=j
                        i=WhichExpansionBasis(1,indice)
                        if(i>0) Hd_HdEC(istate,jstate).Array(i)=ReferencePoint.dH(j,istate,jstate)
                    end do
                end do
            end do
        call InitializeHdLeastSquareFit()
    case('ContinueFitting')
        if(SameTrainingSet) then!Read the rearranged training set
            open(unit=99,file='ReferencePoint.CheckPoint',status='old')
                allocate(ReferencePoint.geom(InternalDimension))
                allocate(ReferencePoint.energy(NState))
                allocate(ReferencePoint.H(NState,NState))
                allocate(ReferencePoint.dH(InternalDimension,NState,NState))
                read(99,*)ReferencePoint.geom
                read(99,*)ReferencePoint.energy
                read(99,*)ReferencePoint.H
                read(99,*)ReferencePoint.dH
            close(99)
            CharTemp128='DegeneratePoint.CheckPoint'
            call ReadDegenerateData(CharTemp128,DegeneratePoint,NDegeneratePoints,InternalDimension)
            NPoints=NPoints-NDegeneratePoints
            CharTemp128='point.CheckPoint'
                allocate(point(NPoints))
                do i=1,NPoints
                    allocate(point(i).geom(InternalDimension))
                    allocate(point(i).energy(NState))
                    allocate(point(i).dH(InternalDimension,NState,NState))
                end do
            call ReadData(CharTemp128,point,NPoints)
            CharTemp128='ArtifactPoint.CheckPoint'
                allocate(ArtifactPoint(NArtifactPoints))
                do i=1,NArtifactPoints
                    allocate(ArtifactPoint(i).geom(InternalDimension))
                    allocate(ArtifactPoint(i).energy(NState))
                end do
            call ReadArtifactData(CharTemp128,ArtifactPoint,NArtifactPoints)
            if(GradualFit) then
                !Read the Cartesian distances to the reference geometry
                allocate(GeomDifference(NPoints))
                open(unit=99,file='GeomDifference.CheckPoint',status='old')
                    do i=1,NPoints
                        read(99,*)GeomDifference(i)
                    end do
                close(99)
            end if
            call InitializeDiabaticHamiltonian(NState,InternalDimension)
        else!Read training set then rearrange it, and check whether reference point has changed
            if(IndexReference==0) then!Same old reference point
                call Initialize_NewTrainingSet()!Read training set then rearrange it
                call InitializeDiabaticHamiltonian(NState,InternalDimension)
            else!Reference point is possibly changed
                !Read old reference point
                allocate(OldRefGeom(InternalDimension)); allocate(OldRefEnergy(NState))
                open(unit=99,file='ReferencePoint.CheckPoint',status='old')
                    read(99,*)OldRefGeom; read(99,*)OldRefEnergy
                close(99)
                call Initialize_NewTrainingSet()!Read training set then rearrange it
                call InitializeDiabaticHamiltonian(NState,InternalDimension)
                flag=.false.!Check whether the reference point has been changed
                do i=1,InternalDimension
                    dbletemp=Abs(ReferencePoint.geom(i)-OldRefGeom(i))
                    if(dbletemp>1d-10.and.dbletemp/Abs(OldRefGeom(i))>1d-10) then
                        flag=.true.; exit
                    end if
                end do
                if(flag) then!Shift Hd expansion origin to new reference point
                    call OriginShift(ReferencePoint.geom-OldRefGeom)
                    i=WhichExpansionBasis(0,indice(1:0))
                    if(i>0) then!Const term shift takes care of the energy zero point shift
                        dbletemp=OldRefEnergy(1)-ReferencePoint.energy(1)
                        forall(istate=1:NState)
                            Hd_HdEC(istate,istate).Array(i)=Hd_HdEC(istate,istate).Array(i)+dbletemp
                        end forall
                    end if
                end if
                deallocate(OldRefGeom); deallocate(OldRefEnergy)!Clean up
            end if
        end if
        call InitializeHdLeastSquareFit()
	case default
		open(unit=99,file='ReferencePoint.CheckPoint',status='old')
            allocate(ReferencePoint.geom(InternalDimension))
            allocate(ReferencePoint.energy(NState))
            allocate(ReferencePoint.H(NState,NState))
            allocate(ReferencePoint.dH(InternalDimension,NState,NState))
            read(99,*)ReferencePoint.geom
            read(99,*)ReferencePoint.energy
            read(99,*)ReferencePoint.H
            read(99,*)ReferencePoint.dH
        close(99)
        call InitializeDiabaticHamiltonian(NState,InternalDimension)
    end select
end subroutine Initialize
subroutine Initialize_NewTrainingSet()!Support Initialize
    character*128::CharTemp128; logical::degenerate
    integer::index,index2,i,istate,jstate,ip
    integer,allocatable,dimension(:)::indices
    real*8,allocatable,dimension(:)::differencetemp
    real*8,dimension(NState)::eigval; real*8,dimension(NState,NState)::eigvec
	real*8,dimension(CartesianDimension)::grad1,grad2,g,h
    type(Data)::ReferencePointtemp
    type(Data),allocatable,dimension(:)::pointtemp,pointswap,ArtifactPointtemp
    !Read training set
        !Allocate storage space
            allocate(pointtemp(NPoints))
            do ip=1,NPoints
                allocate(pointtemp(ip).geom(CartesianDimension))
                allocate(pointtemp(ip).energy(NState))
                allocate(pointtemp(ip).dH(CartesianDimension,NState,NState))
            end do
        call ReadElectronicStructureData(pointtemp,NPoints)
        allocate(ArtifactPointtemp(NArtifactPoints))
        open(unit=99,file=ArtifactGeometryDataFile,status='old',iostat=istate)
		open(unit=100,file=ArtifactEnergyDataFile,status='old',iostat=jstate)
		    if(istate==0.and.jstate==0) then
                do ip=1,NArtifactPoints
                    allocate(ArtifactPointtemp(ip).geom(CartesianDimension))
                    read(99,*)ArtifactPointtemp(ip).geom
                    allocate(ArtifactPointtemp(ip).energy(NState))
                    read(100,*)ArtifactPointtemp(ip).energy
				end do
			else
				if(NArtifactPoints/=0) stop 'Program abort: artifact data not found'
			end if
        close(100)
        close(99)
    !Obtain the reference point
        !Allocate storage space
            allocate(ReferencePointtemp.geom(CartesianDimension))
            allocate(ReferencePointtemp.energy(NState))
            allocate(ReferencePointtemp.dH(CartesianDimension,NState,NState))
            allocate(ReferencePoint.geom(InternalDimension))
            allocate(ReferencePoint.energy(NState))
            allocate(ReferencePoint.H(NState,NState))
            allocate(ReferencePoint.dH(InternalDimension,NState,NState))
        if(IndexReference==0) then!Use old reference point
            open(unit=99,file='ReferencePoint.CheckPoint',status='old')!Obain in internal coordinate
                read(99,*)ReferencePoint.geom
                read(99,*)ReferencePoint.energy
                read(99,*)ReferencePoint.H
                read(99,*)ReferencePoint.dH
            close(99)
            !Convert reference point from internal coordinate to Cartesian coordinate
            call Internal2Cartesian(ReferencePoint.geom,ReferencePoint.dH,ReferencePointtemp.geom,ReferencePointtemp.dH,InternalDimension,CartesianDimension,NState,r0=reshape(MoleculeDetail.RefConfig,[CartesianDimension]))
        else!Use new reference point
            ReferencePointtemp=pointtemp(IndexReference)!Obtain in Cartesian coordinate
            !Convert reference point from Cartesian coordinate to internal coordinate, and transform if degenerate
            ReferencePoint.energy=ReferencePointtemp.energy
            call Cartesian2Internal(ReferencePointtemp.geom,ReferencePointtemp.dH,ReferencePoint.geom,ReferencePoint.dH,CartesianDimension,InternalDimension,NState)
            call CheckDegeneracy(degenerate,AlmostDegenerate,ReferencePoint.energy,NState)
            if(Degenerate) then
                call NondegenerateRepresentation(ReferencePoint.dH,eigval,eigvec,InternalDimension,NState,DegenerateThreshold=AlmostDegenerate)
                ReferencePoint.H=transpose(eigvec)
                forall(istate=1:NState)
                    ReferencePoint.H(:,istate)=ReferencePoint.energy(istate)*ReferencePoint.H(:,istate)
                end forall
                ReferencePoint.H=matmul(ReferencePoint.H,eigvec)-ReferencePoint.energy(1)*UnitMatrix(NState)
            end if
        end if
    !Rearrange training set: shift energy zero point to the ground state energy of the reference point
    !                        scale fitting weight according to the ground state energy of the point
    do ip=1,NPoints!Modify points
        pointtemp(ip).energy=pointtemp(ip).energy-ReferencePoint.energy(1)
        if(pointtemp(ip).energy(1)>HighEnergy) pointtemp(ip).weight=(HighEnergy/pointtemp(ip).energy(1))**2
    end do
    !If want to perform the fitting procedure gradually, we will require:
    !    All geometries aligned
    !    The Cartesian distances to the reference geometry
    !    All data points sorted according to the Cartesian distance to the reference geometry in ascending order
    if(GradualFit) then
        allocate(GeomDifference(NPoints))
        if(Unaligned) then!Align all geometries and calculate Cartesian distances to the reference geometry
            call StandardizeGeometry(ReferencePointtemp.geom,MoleculeDetail.mass,MoleculeDetail.NAtoms,NState,&
                grad=ReferencePointtemp.dH)
            do ip=1,NPoints
                call StandardizeGeometry(pointtemp(ip).geom,MoleculeDetail.mass,MoleculeDetail.NAtoms,NState,&
                ref=ReferencePointtemp.geom,diff=GeomDifference(ip),grad=pointtemp(ip).dH)
            end do
        else!Directly calculate the Cartesian distances
            do ip=1,NPoints
                GeomDifference(ip)=dot_product(pointtemp(ip).geom-ReferencePointtemp.geom,pointtemp(ip).geom-ReferencePointtemp.geom)
            end do
        end if
        if(Unsorted) then!Sort all data points according to the Cartesian distance ascendingly
            !Allocate work space
                allocate(indices(NPoints)); forall(i=1:NPoints); indices(i)=i; end forall
                allocate(pointswap(NPoints))
                do ip=1,NPoints
                    allocate(pointswap(ip).geom(CartesianDimension))
                    allocate(pointswap(ip).energy(NState))
                    allocate(pointswap(ip).dH(CartesianDimension,NState,NState))
                end do
            call dQuickSort(GeomDifference,1,NPoints,indices,NPoints)
            do ip=1,NPoints; pointswap(ip)=pointtemp(indices(ip)); end do; pointtemp=pointswap
            deallocate(indices)!Clean up
            do ip=1,NPoints
                deallocate(pointswap(ip).geom); deallocate(pointswap(ip).energy); deallocate(pointswap(ip).dH)
            end do; deallocate(pointswap)
        end if
        open(unit=99,file='GeomDifference.CheckPoint',status='replace')!Save the Cartesian distances
            do ip=1,NPoints; write(99,*)GeomDifference(ip); end do
        close(99)
    else if(Unaligned) then!Also align all geometries for future Cartesian coordinate analysis
        call StandardizeGeometry(ReferencePointtemp.geom,MoleculeDetail.mass,MoleculeDetail.NAtoms,NState,&
            grad=ReferencePointtemp.dH)
        do ip=1,NPoints
            call StandardizeGeometry(pointtemp(ip).geom,MoleculeDetail.mass,MoleculeDetail.NAtoms,NState,&
            ref=ReferencePointtemp.geom,grad=pointtemp(ip).dH)
        end do
    end if
	!Provide a human readable version of training set
		open(unit=99,file='TrainingEnergy.txt',status='replace')
		    write(99,'(A10)',advance='no')'Geometry#'//char(9)
		    do istate=1,NState-1
		    	write(99,'(2x,A6,I2,A5,3x,A1)',advance='no')'Energy',istate,'/cm-1',char(9)
		    end do
		    write(99,'(2x,A6,I2,A5,3x)')'Energy',NState,'/cm-1'
		    do ip=1,NPoints
		    	write(99,'(I9,A1)',advance='no')ip,char(9)
		    	do jstate=1,NState-1
		    		write(99,'(F18.8,A1)',advance='no')pointtemp(ip).energy(jstate)/cm_1InAU,char(9)
		    	end do
		    	write(99,'(F18.8)')pointtemp(ip).energy(NState)/cm_1InAU
		    end do
		close(99)
		open(unit=99,file='TrainingGradient.txt',status='replace')
            write(99,'(A10)',advance='no')'Geometry#'//char(9)
            do istate=1,NState!Energy gradient label
                write(99,'(A11,I2,A5,A1)',advance='no')'Energy grad',istate,'/a.u.',char(9)
			end do
			do istate=1,NState!Interstate coupling label
        		do jstate=istate+1,NState
        			write(99,'(2x,A3,I2,A1,I2,A5,3x,A1)',advance='no')'ISC',jstate,'&',istate,'/a.u.',char(9)
        		end do
            end do
        	do istate=1,NState!Nonadiabatic coupling label
        		do jstate=istate+1,NState
        			write(99,'(2x,A3,I2,A1,I2,A5,3x,A1)',advance='no')'NAC',jstate,'&',istate,'/a.u.',char(9)
        		end do
            end do
            write(99,*)
            do ip=1,NPoints
                call CheckDegeneracy(degenerate,1d-10,pointtemp(ip).energy,NState)
			    if(degenerate) then!gh orthogonalization for conical intersection
			    	do istate=1,NState-1!Identify the degenerate states
			    		if(pointtemp(ip).energy(istate+1)-pointtemp(ip).energy(istate)<1d-10) exit
			    	end do
			    	grad1=pointtemp(ip).dH(:,istate,istate)
			    	grad2=pointtemp(ip).dH(:,istate+1,istate+1)
			    	h    =pointtemp(ip).dH(:,istate+1,istate)
			    	call ghOrthogonalization(grad1,grad2,h,CartesianDimension)
			    	g=(grad2-grad1)/2d0
			    	write(CharTemp128,*)ip
			    	open(unit=100,file='g'//trim(adjustl(CharTemp128)),status='replace')
			    	    write(100,*)g
			    	close(100)
			    	open(unit=100,file='h'//trim(adjustl(CharTemp128)),status='replace')
			    	    write(100,*)h
			    	close(100)
			    end if
                write(99,'(I9,A1)',advance='no')ip,char(9)
                do jstate=1,NState!Energy gradient
                    write(99,'(F18.8,A1)',advance='no')norm2(pointtemp(ip).dH(:,jstate,jstate)),char(9)
				end do
				do istate=1,NState!Interstate coupling
        			do jstate=istate+1,NState
        				write(99,'(F18.8,A1)',advance='no')norm2(pointtemp(ip).dH(:,jstate,istate)),char(9)
        			end do
        		end do
        		do istate=1,NState!Nonadiabatic coupling
        			do jstate=istate+1,NState
        				write(99,'(F18.8,A1)',advance='no')norm2(pointtemp(ip).dH(:,jstate,istate))/dABS(pointtemp(ip).energy(istate)-pointtemp(ip).energy(jstate)),char(9)
        			end do
        		end do
                write(99,*)
            end do
		close(99)
    do ip=1,NArtifactPoints!Modify artifact points
        ArtifactPointtemp(ip).energy=ArtifactPointtemp(ip).energy-ReferencePoint.energy(1)
        if(ArtifactPointtemp(ip).energy(1)>HighEnergy) ArtifactPointtemp(ip).weight=(HighEnergy/ArtifactPointtemp(ip).energy(1))**2
    end do
    !Convert points from Cartesian coordinate to internal coordinate, and identify almost degenerate points
        allocate(pointswap(NPoints))
        do ip=1,NPoints
            pointswap(ip).weight=pointtemp(ip).weight
            allocate(pointswap(ip).geom(InternalDimension))
            allocate(pointswap(ip).energy(NState))
                pointswap(ip).energy=pointtemp(ip).energy
            allocate(pointswap(ip).dH(InternalDimension,NState,NState))
            call Cartesian2Internal(pointtemp(ip).geom,pointtemp(ip).dH,pointswap(ip).geom,pointswap(ip).dH,CartesianDimension,InternalDimension,NState)
            pointswap(ip).geom=pointswap(ip).geom-ReferencePoint.geom!This program requires only internal coordinate difference
        end do
        call IdentifyDegeneracy(DegeneratePoint,NDegeneratePoints,indices,pointswap,NPoints)
        NPoints=NPoints-NDegeneratePoints
        allocate(point(NPoints))
        index=1
        index2=1
        do ip=1,NPoints+NDegeneratePoints
            if(ip==indices(index)) then
                index=index+1
            else
                allocate(point(index2).geom(InternalDimension))
                allocate(point(index2).energy(NState))
                allocate(point(index2).dH(InternalDimension,NState,NState))
                point(index2)=pointswap(ip)
                index2=index2+1
            end if
        end do
        if(GradualFit) then!Almost degenerate points are kicked out from gradual fit procedure
            allocate(differencetemp(NPoints+NDegeneratePoints))
            differencetemp=GeomDifference
            deallocate(GeomDifference)
            allocate(GeomDifference(NPoints))
            index=1
            index2=1
            do ip=1,NPoints+NDegeneratePoints
                if(ip==indices(index)) then
                    index=index+1
                else
                    GeomDifference(index2)=differencetemp(ip)
                    index2=index2+1
                end if
            end do
            !Write down the Cartesian distances to the reference geometry
            open(unit=99,file='GeomDifference.CheckPoint',status='replace')
                do ip=1,NPoints
                    write(99,*)GeomDifference(ip)
                end do
            close(99)
        end if
    !Convert artifact points from Cartesian coordinate to internal coordinate
        allocate(ArtifactPoint(NArtifactPoints))
        do ip=1,NArtifactPoints
            ArtifactPoint(ip).weight=ArtifactPointtemp(ip).weight
            allocate(ArtifactPoint(ip).geom(InternalDimension))
            allocate(ArtifactPoint(ip).energy(NState))
                ArtifactPoint(ip).energy=ArtifactPointtemp(ip).energy
            !Artifact points do not have energy gradient and interstate coupling
            call InternalCoordinate(ArtifactPointtemp(ip).geom,ArtifactPoint(ip).geom,CartesianDimension,InternalDimension)
            ArtifactPoint(ip).geom=ArtifactPoint(ip).geom-ReferencePoint.geom!This program requires only internal coordinate difference
        end do
    !Clean up
        deallocate(indices)
        deallocate(ReferencePointtemp.geom); deallocate(ReferencePointtemp.energy); deallocate(ReferencePointtemp.dH)
        do ip=1,NPoints
            deallocate(pointtemp(ip).geom); deallocate(pointtemp(ip).energy); deallocate(pointtemp(ip).dH)
        end do; deallocate(pointtemp)
        do ip=1,NPoints
            deallocate(pointswap(ip).geom); deallocate(pointswap(ip).energy); deallocate(pointswap(ip).dH)
        end do; deallocate(pointswap)
        do ip=1,NArtifactPoints
            deallocate(ArtifactPointtemp(ip).geom); deallocate(ArtifactPointtemp(ip).energy)
        end do; deallocate(ArtifactPointtemp)
    !Save the rearranged training set
    open(unit=99,file='ReferencePoint.CheckPoint',status='replace')
        write(99,*)ReferencePoint.geom
        write(99,*)ReferencePoint.energy
        write(99,*)ReferencePoint.H
        write(99,*)ReferencePoint.dH
    close(99)
    CharTemp128='point.CheckPoint'
    call WriteData(CharTemp128,point,NPoints)
    CharTemp128='DegeneratePoint.CheckPoint'
    call WriteDegenerateData(CharTemp128,DegeneratePoint,NDegeneratePoints)
    CharTemp128='ArtifactPoint.CheckPoint'
    call WriteArtifactData(CharTemp128,ArtifactPoint,NArtifactPoints)
end subroutine Initialize_NewTrainingSet

end program main