!SurfGenBound: surface generation package for bounded system
!Construct diabatic hamiltonian (Hd) by least square fitting H & ▽H
!IO & computation unit: atomic unit
!
!Level: Main <- (HdLeastSquareFit, Analyze, NadVibS) <- ElectronicStructure <- Basic <- DiabaticHamiltonian
program main
    use Basic
    use ElectronicStructure
    use HdLeastSquareFit
    use Analyzation
    use NadVibS
    implicit none
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
    write(*,*)'Electronic structure software in use = '//ElectronicStructureSoftware
    call ReadInput()
    call Initialize()
!------------- End --------------

!----------- Run job ------------
    call showtime()
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
                NPoints=ceiling(dble(NExpansionCoefficients)/dble(DataPerPoint))
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
        case('Analyze')
            call Analyze()
        case('NadVibS')
            call GenerateNadVibSInput()
        case default!Throw a warning
            write(*,'(1x,A35,1x,A32)')'Program abort: unsupported job type',JobType
            stop
    end select
!------------- End --------------

!---------- Clean up ------------
    write(*,*)'Mission complete'
!------------- End --------------

contains
    !Read main input files: input, MoleculeDetailFile, AdvancedInput (optional)
    subroutine ReadInput()
        character*128::MoleculeDetailFile
        logical::advance
        integer::i
        !Read main input, write some job comment
        open(unit=99,file='input',status='old')
            read(99,*)
            read(99,*)
            read(99,*)
            read(99,*)JobType
                write(*,*)'Job type: '//JobType
            read(99,*)
            read(99,*)MoleculeDetailFile
            read(99,*)
            read(99,*)NStates
            read(99,*)
            read(99,*)NOrder
            read(99,*)
            read(99,*)SameTrainingSet
            read(99,*)
            read(99,*)NPoints
            read(99,*)
            read(99,*)IndexReference
            read(99,*)
            read(99,*)ArtifactGeometryDataFile
            read(99,*)
            read(99,*)ArtifactEnergyDataFile
            read(99,*)
            read(99,*)NArtifactPoints
            read(99,*)
            read(99,*)advance
        close(99)
        open(unit=99,file=MoleculeDetailFile,status='old')
            read(99,*)
            read(99,*)
            read(99,*)
            read(99,*)NAtoms
                allocate(MoleculeDetail.ElementSymbol(NAtoms))
                allocate(MoleculeDetail.RefConfig(3,NAtoms))
                allocate(MoleculeDetail.mass(NAtoms))
            read(99,*)
            do i=1,NAtoms
                read(99,'(A2,3F20.15)')MoleculeDetail.ElementSymbol(i),MoleculeDetail.RefConfig(:,i)
            end do
            read(99,*)
            do i=1,NAtoms
                read(99,*)MoleculeDetail.mass(i)
            end do
                MoleculeDetail.mass=MoleculeDetail.mass*AMUInAU!Convert to atomic unit
        close(99)
        if(advance) then!If requested, read advanced input
            write(*,*)'Advanced input requested, parameters are set to user specification'
            open(unit=99,file='advance.in',status='old')
                namelist /AdvancedInput/ &
                    !Basic
                        HighEnergy,AlmostDegenerate,&
                    !HdLeastSquareFit
                        HdLSF_Regularization,HdLSF_Solver,&
                        HdLSF_MaxHopperIteration,HdLSF_MaxLocalMinimizerIteration,HdLSF_Max2StepIteration,&
                        HdLSF_pseudolinearFollowFreq,HdLSF_pseudolinearMaxMonotonicalIncrease,&
                        HdLSF_LineSearcher,HdLSF_UseStrongWolfe,HdLSF_LBFGSMemory,HdLSF_ConjugateGradientSolver,&
                    !Main
                        GradualFit,AutoGradualFit,ManualNPoints,Unaligned,Unsorted
                read(99,nml=AdvancedInput)
            close(99)
        end if
    end subroutine ReadInput

    !The initializer for the program
    subroutine Initialize()
        character*128::CharTemp128
        integer::ip
        real*8::absdev
        real*8,allocatable,dimension(:)::OldRefGeom
        !General initialize
            call BetterRandomSeed()
            call InitializeBasic()
            CartesianDimension=3*NAtoms
            call DefineInternalCoordinate(ElectronicStructureSoftware,InternalDimension)
        select case(JobType)!Job specific initialize
            case('FitNewDiabaticHamiltonian')!To fit Hd from scratch, read training set then rearrange it
                call Initialize_NewTrainingSet()
                call InitializeDiabaticHamiltonian()!Must call after reference point is known
                call InitializeHdLeastSquareFit()
            case('ContinueFitting')
                if(SameTrainingSet) then!Read the rearranged training set
                    open(unit=99,file='ReferencePoint.CheckPoint',status='old')
                        allocate(ReferencePoint.geom(InternalDimension))
                        allocate(ReferencePoint.energy(NStates))
                        allocate(ReferencePoint.H(NStates,NStates))
                        allocate(ReferencePoint.dH(InternalDimension,NStates,NStates))
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
                        do ip=1,NPoints
                            allocate(point(ip).geom(InternalDimension))
                            allocate(point(ip).energy(NStates))
                            allocate(point(ip).dH(InternalDimension,NStates,NStates))
                        end do
                    call ReadData(CharTemp128,point,NPoints)
                    CharTemp128='ArtifactPoint.CheckPoint'
                        allocate(ArtifactPoint(NArtifactPoints))
                        do ip=1,NArtifactPoints
                            allocate(ArtifactPoint(ip).geom(InternalDimension))
                            allocate(ArtifactPoint(ip).energy(NStates))
                        end do
                    call ReadArtifactData(CharTemp128,ArtifactPoint,NArtifactPoints)
                    if(GradualFit) then
                        !Read the Cartesian distances to the reference geometry
                        allocate(GeomDifference(NPoints))
                        open(unit=99,file='GeomDifference.CheckPoint',status='old')
                            do ip=1,NPoints
                                read(99,*)GeomDifference(ip)
                            end do
                        close(99)
                    end if
                else!Read training set then rearrange it, and check whether reference point has changed
                    call Initialize_NewTrainingSet()
                    !Check whether the reference point has been changed
                    allocate(OldRefGeom(InternalDimension))
                    open(unit=99,file='ReferencePoint.CheckPoint',status='old')
                        read(99,*)OldRefGeom
                    close(99)
                    do ip=1,InternalDimension
                        absdev=Abs(ReferencePoint.geom(ip)-OldRefGeom(ip))
                        if(absdev>1d-14.and.absdev/Abs(OldRefGeom(ip))>1d-14) then
                            ReferenceChange=.true.
                            exit
                        end if
                    end do
                end if
                call InitializeDiabaticHamiltonian()
                call InitializeHdLeastSquareFit()
            case default
                call InitializeDiabaticHamiltonian()
        end select
    end subroutine Initialize
    !Support Initialize
    subroutine Initialize_NewTrainingSet()
        character*128::CharTemp128
        logical::degenerate
        integer::index,index2,i,istate,jstate,ip
        integer,allocatable,dimension(:)::indices
        real*8,allocatable,dimension(:)::differencetemp
        real*8,dimension(NStates)::eigval
        real*8,dimension(NStates,NStates)::eigvec
        type(Data)::ReferencePointtemp
        type(Data),allocatable,dimension(:)::pointtemp,pointswap,ArtifactPointtemp
        !Read training set
            !Allocate storage space
                allocate(pointtemp(NPoints))
                do ip=1,NPoints
                    allocate(pointtemp(ip).geom(CartesianDimension))
                    allocate(pointtemp(ip).energy(NStates))
                    allocate(pointtemp(ip).dH(CartesianDimension,NStates,NStates))
                end do
            call ReadElectronicStructureData(pointtemp,NPoints)
            allocate(ArtifactPointtemp(NArtifactPoints))
            open(unit=99,file=ArtifactGeometryDataFile,status='old')
            open(unit=100,file=ArtifactEnergyDataFile,status='old')
                do ip=1,NArtifactPoints
                    allocate(ArtifactPointtemp(ip).geom(CartesianDimension))
                    read(99,*)ArtifactPointtemp(ip).geom
                    allocate(ArtifactPointtemp(ip).energy(NStates))
                    read(100,*)ArtifactPointtemp(ip).energy
                end do
            close(100)
            close(99)
        !Rearrange training set: shift energy zero point to the ground state energy of the reference point
        !                        scale fitting weight according to the ground state energy of the point
        !The reference point
            !Allocate storage space
                allocate(ReferencePointtemp.geom(CartesianDimension))
                allocate(ReferencePointtemp.energy(NStates))
                allocate(ReferencePointtemp.dH(CartesianDimension,NStates,NStates))
            ReferencePointtemp=pointtemp(IndexReference)
        do ip=1,NPoints!Modify points
            pointtemp(ip).energy=pointtemp(ip).energy-ReferencePointtemp.energy(1)
            if(pointtemp(ip).energy(1)>HighEnergy) pointtemp(ip).weight=HighEnergy/pointtemp(ip).energy(1)
        end do
        do ip=1,NArtifactPoints!Modify artifact points
            ArtifactPointtemp(ip).energy=ArtifactPointtemp(ip).energy-ReferencePointtemp.energy(1)
            if(ArtifactPointtemp(ip).energy(1)>HighEnergy) ArtifactPointtemp(ip).weight=HighEnergy/ArtifactPointtemp(ip).energy(1)
        end do
        !If want to perform the fitting procedure gradually, we will require:
        !    All geometries aligned
        !    The Cartesian distances to the reference geometry
        !    All data points sorted according to the Cartesian distance to the reference geometry in ascending order
        if(GradualFit) then
            allocate(GeomDifference(NPoints))
            if(Unaligned) then
                call StandardizeGeometry(ReferencePointtemp.geom,MoleculeDetail.mass,NAtoms,NStates,&
                    nadgrad=ReferencePointtemp.dH)
                do ip=1,NPoints
                    call StandardizeGeometry(pointtemp(ip).geom,MoleculeDetail.mass,NAtoms,NStates,&
                        nadgrad=pointtemp(ip).dH,reference=ReferencePointtemp.geom,difference=GeomDifference(ip))
                end do
            else
                do ip=1,NPoints
                    GeomDifference(ip)=dot_product(pointtemp(ip).geom-ReferencePointtemp.geom,pointtemp(ip).geom-ReferencePointtemp.geom)
                end do
            end if
            if(Unsorted) then
                !Allocate work space
                    allocate(indices(NPoints))
                    forall(i=1:NPoints)
                        indices(i)=i
                    end forall
                    allocate(pointswap(NPoints))
                    do ip=1,NPoints
                        allocate(pointswap(ip).geom(CartesianDimension))
                        allocate(pointswap(ip).energy(NStates))
                        allocate(pointswap(ip).dH(CartesianDimension,NStates,NStates))
                    end do
                call dQuickSort(GeomDifference,1,NPoints,indices,NPoints)
                do ip=1,NPoints
                    pointswap(ip)=pointtemp(indices(ip))
                end do
                pointtemp=pointswap
                !Clean up
                    deallocate(indices)
                    deallocate(pointswap)
            end if
            !Write down the Cartesian distances to the reference geometry
            open(unit=99,file='GeomDifference.CheckPoint',status='replace')
                do ip=1,NPoints
                    write(99,*)GeomDifference(ip)
                end do
            close(99)
        end if
        !This program requires only internal coordinate difference
        !Convert reference point from Cartesian coordinate to internal coordinate, and transform if degenerate
            allocate(ReferencePoint.geom(InternalDimension))
            allocate(ReferencePoint.energy(NStates))
                ReferencePoint.energy=ReferencePointtemp.energy
            allocate(ReferencePoint.H(NStates,NStates))
            allocate(ReferencePoint.dH(InternalDimension,NStates,NStates))
            call Cartesian2Internal(ReferencePointtemp.geom,CartesianDimension,ReferencePoint.geom,InternalDimension,NStates,&
                cartnadgrad=ReferencePointtemp.dH,intnadgrad=ReferencePoint.dH)
            call CheckDegeneracy(degenerate,AlmostDegenerate,ReferencePoint.energy,NStates)
            if(Degenerate) then
                call NondegenerateRepresentation(ReferencePoint.dH,eigval,eigvec,InternalDimension,NStates)
                ReferencePoint.H=transpose(eigvec)
                forall(istate=1:NStates)
                    ReferencePoint.H(:,istate)=ReferencePoint.energy(istate)*ReferencePoint.H(:,istate)
                end forall
                ReferencePoint.H=matmul(ReferencePoint.H,eigvec)-ReferencePoint.energy(1)*UnitMatrix(NStates)
            end if
        !Convert points from Cartesian coordinate to internal coordinate, and identify almost degenerate points
            allocate(pointswap(NPoints))
            do ip=1,NPoints
                pointswap(ip).weight=pointtemp(ip).weight
                allocate(pointswap(ip).geom(InternalDimension))
                allocate(pointswap(ip).energy(NStates))
                    pointswap(ip).energy=pointtemp(ip).energy
                allocate(pointswap(ip).dH(InternalDimension,NStates,NStates))
                call Cartesian2Internal(pointtemp(ip).geom,CartesianDimension,pointswap(ip).geom,InternalDimension,NStates,&
                    cartnadgrad=pointtemp(ip).dH,intnadgrad=pointswap(ip).dH)
                pointswap(ip).geom=pointswap(ip).geom-ReferencePoint.geom
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
                    allocate(point(index2).energy(NStates))
                    allocate(point(index2).dH(InternalDimension,NStates,NStates))
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
                allocate(ArtifactPoint(ip).energy(NStates))
                    ArtifactPoint(ip).energy=ArtifactPointtemp(ip).energy
                !Artifact points do not have energy gradient and interstate coupling
                ArtifactPoint(ip).geom=InternalCoordinateq(ArtifactPoint(ip).geom,InternalDimension,CartesianDimension)-ReferencePoint.geom
            end do
        !Clean up
            deallocate(indices)
            deallocate(pointtemp)
            deallocate(pointswap)
            deallocate(ArtifactPointtemp)
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