!Library interface
!Global parameters, derived types, data storage
!Basic routines
module Basic
    use General
    use Mathematics
    use LinearAlgebra
    use GeometryTransformation
    use NonlinearOptimization
    use Nonadiabatic
    implicit none

!Parameter
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
    !Points with ground state energy difference to the reference geometry higher than this threshold will have lower weight:
    !    weight = HighEnergy / ground state energy difference to the reference geometry
    real*8::HighEnergy=1d2
    !States with energy difference below AlmostDegenerate will be considered degenerate
    real*8::AlmostDegenerate=1d-4

!Derived type
    !Example: type(Data),allocatable,dimension(:)::point, point(ip) stands for the ip-th data point
    type Data
        real*8::weight=1d0!Fitting weight
        real*8,allocatable,dimension(:)::geom,energy
        real*8,allocatable,dimension(:,:)::H
        real*8,allocatable,dimension(:,:,:)::dH!dH(m,i,j) = ▽_m H_i,j
    end type Data

    !See examples of my special .xyz for the meaning of each variable
    !Store the details of an input molecule
    type MoleculeDetails
        character*2,allocatable,dimension(:)::ElementSymbol!Element symbol of each atom
        real*8,allocatable,dimension(:,:)::RefConfig!Short for REFerence CONFIGuration
        real*8,allocatable,dimension(:)::mass!Mass of each atom
    end type MoleculeDetails
    
!Input variable
    !See input example for the meaning of each variable
    !Main input
        character*32::JobType
        integer::NStates,NOrder
        logical::SameTrainingSet
        integer::NPoints,IndexReference
        character*128::ArtifactGeometryDataFile,ArtifactEnergyDataFile
        integer::NArtifactPoints
    !Molecule detail input
        integer::NAtoms
        type(MoleculeDetails)::MoleculeDetail

!Global variable
    integer::NDegeneratePoints!Number of almost degenerate points
    integer::CartesianDimension,InternalDimension!Dimension of Cartesian/internal coordinate space
    !Training set:
        logical::ReferenceChange=.false.
        integer::DataPerPoint,DataPerDegeneratePoint!How much data a point provides
        type(Data)::ReferencePoint!The reference point
        !Data to be least square fitted:
        type(Data),allocatable,dimension(:)::point,&!Regular data point, in adiabatic representation
            DegeneratePoint,&!In other nondegenerate representation
            ArtifactPoint!Unreliable data, do not fit ▽H
    real*8,allocatable,dimension(:)::GeomDifference!Difference between data points geometries and reference point geometry

contains
!The initializer for Basic module
subroutine InitializeBasic()
    call InitializePhaseFixing(NStates)
end subroutine InitializeBasic

!Identify the almost degenerate points, then store them in degpoint
!Input:  unallocated array degpoint, NPoints dimensional array point
!Output:   degpoint harvests almost degenerate data points
!        NDegpoints harvests the  number of almost degenerate data points
!        IndicesDeg harvests the indices of almost degenerate data points in array point
subroutine IdentifyDegeneracy(degpoint,NDegpoints,IndicesDeg,point,NPoints)
    type(Data),allocatable,dimension(:),intent(out)::degpoint
    integer,intent(out)::NDegpoints
    integer,allocatable,dimension(:),intent(out)::IndicesDeg
    integer,intent(in)::NPoints
    type(Data),dimension(NPoints),intent(in)::point
    logical::degenerate
    integer::dim,index,ip,istate,jstate
    real*8,dimension(NStates)::eigval
    real*8,dimension(NStates,NStates)::eigvec
    !Count how many points are almost degenerate
    NDegpoints=0
    do ip=1,NPoints
        call CheckDegeneracy(degenerate,AlmostDegenerate,point(ip).energy,NStates)
        if(degenerate) NDegpoints=NDegpoints+1
    end do
    !Tag the almost degenerate points, and transform them into nondegenerate representation
    if(allocated(IndicesDeg)) deallocate(IndicesDeg)
    allocate(IndicesDeg(NDegpoints))
    if(allocated(degpoint)) deallocate(degpoint)
    allocate(degpoint(NDegpoints))
    dim=size(point(1).geom)
    index=1
    do ip=1,NPoints
        call CheckDegeneracy(degenerate,AlmostDegenerate,point(ip).energy,NStates)
        if(degenerate) then
            IndicesDeg(index)=ip
            degpoint(index).weight=point(ip).weight
            allocate(degpoint(index).geom(dim))
            degpoint(index).geom=point(ip).geom
            allocate(degpoint(index).energy(dim))
            degpoint(index).energy=point(ip).energy
            allocate(degpoint(index).dH(dim,NStates,NStates))
            degpoint(index).dH=point(ip).dH
            call NondegenerateRepresentation(degpoint(index).dH,eigval,eigvec,dim,NStates)
            allocate(degpoint(index).H(NStates,NStates))
            degpoint(index).H=transpose(eigvec)
            forall(istate=1:NStates)
                degpoint(index).H(:,istate)=degpoint(index).energy(istate)*degpoint(index).H(:,istate)
            end forall
            degpoint(index).H=matmul(degpoint(index).H,eigvec)
            index=index+1
        end if
    end do
end subroutine IdentifyDegeneracy

!-------------------- Input / Output --------------------
    !Load data in DataFile to point
    subroutine ReadData(DataFile,point,NPoints)
        character*128,intent(in)::DataFile
        integer,intent(in)::NPoints
        type(Data),dimension(NPoints),intent(inout)::point
        integer::ip,istate,jstate
        open(unit=99,file=DataFile,status='old')
            do ip=1,NPoints
                read(99,*)point(ip).weight
                read(99,*)point(ip).geom
                read(99,*)point(ip).energy
                do istate=1,NStates
                    do jstate=istate,NStates
                        read(99,*)point(ip).dH(:,jstate,istate)
                    end do
                end do
            end do
        close(99)
    end subroutine ReadData
    
    !Write data in point to DataFile
    subroutine WriteData(DataFile,point,NPoints)
        character*128,intent(in)::DataFile
        integer,intent(in)::NPoints
        type(Data),dimension(NPoints),intent(in)::point
        integer::ip,istate,jstate
        open(unit=99,file=DataFile,status='replace')
            do ip=1,NPoints
                write(99,*)point(ip).weight
                write(99,*)point(ip).geom
                write(99,*)point(ip).energy
                do istate=1,NStates
                    do jstate=istate,NStates
                        write(99,*)point(ip).dH(:,jstate,istate)
                    end do
                end do
            end do
        close(99)
    end subroutine WriteData

    !Load degenerate data in DataFile to point
    !Input:  unallocated array point, dim is the dimension of coordinate space
    !Output: point harvests degenerate data in DataFile, NPoints harvests the number of degenerate data points
    subroutine ReadDegenerateData(DataFile,point,NPoints,dim)
        character*128,intent(in)::DataFile
        type(Data),allocatable,dimension(:),intent(out)::point
        integer,intent(out)::NPoints
        integer,intent(in)::dim
        integer::ip,istate,jstate
        open(unit=99,file=DataFile,status='old')
            read(99,*)NPoints
            allocate(point(NPoints))
            do ip=1,NPoints
                read(99,*)point(ip).weight
                allocate(point(ip).geom(dim))
                read(99,*)point(ip).geom
                allocate(point(ip).energy(NStates))
                read(99,*)point(ip).energy
                allocate(point(ip).H(NStates,NStates))
                read(99,*)point(ip).H
                allocate(point(ip).dH(dim,NStates,NStates))
                do istate=1,NStates
                    do jstate=istate,NStates
                        read(99,*)point(ip).dH(:,jstate,istate)
                    end do
                end do
            end do
        close(99)
    end subroutine ReadDegenerateData

    !Write degenerate data in point to DataFile
    subroutine WriteDegenerateData(DataFile,point,NPoints)
        character*128,intent(in)::DataFile
        integer,intent(in)::NPoints
        type(Data),dimension(NPoints),intent(in)::point
        integer::ip,istate,jstate
        open(unit=99,file=DataFile,status='replace')
            write(99,*)NPoints
            do ip=1,NPoints
                write(99,*)point(ip).weight
                write(99,*)point(ip).geom
                write(99,*)point(ip).energy
                write(99,*)point(ip).H
                do istate=1,NStates
                    do jstate=istate,NStates
                        write(99,*)point(ip).dH(:,jstate,istate)
                    end do
                end do
            end do
        close(99)
    end subroutine WriteDegenerateData
    
    !Load artifact data in DataFile to point
    subroutine ReadArtifactData(DataFile,point,NPoints)
        character*128,intent(in)::DataFile
        integer,intent(in)::NPoints
        type(Data),dimension(NPoints),intent(inout)::point
        integer::ip
        open(unit=99,file=DataFile,status='old')
            do ip=1,NPoints
                read(99,*)point(ip).weight
                read(99,*)point(ip).geom
                read(99,*)point(ip).energy
            end do
        close(99)
    end subroutine ReadArtifactData
    
    !Write artifact data in point to DataFile
    subroutine WriteArtifactData(DataFile,point,NPoints)
        character*128,intent(in)::DataFile
        integer,intent(in)::NPoints
        type(Data),dimension(NPoints),intent(in)::point
        integer::ip
        open(unit=99,file=DataFile,status='replace')
            do ip=1,NPoints
                write(99,*)point(ip).weight
                write(99,*)point(ip).geom
                write(99,*)point(ip).energy
            end do
        close(99)
    end subroutine WriteArtifactData
!------------------------- End --------------------------

!---------- Other nondegenerate representation ----------
    !When Hamiltonian is almost degenerate, adiabatic basis experiences numerical difficulty
    !    This causes no trouble to H, but operators other than H suffer a lot (e.g. ▽H)
    !Diabatz is always good, but it cannot be determined without knowledge of neighbourhood
    !So we want some other nondegenerate representation based only on single point information
    !Choose (▽H)^2 representation: it is good enough for fitting a C1 molecule though does not incorporate symmetry
    !    (▽H)^2 is Hermitian so there exists exactly one (▽H)^2 representation
    !    It reduces to s orthogonal to h for 2-fold degeneracy case, where s is the average force of 2 PESs
    !Note it is only appropriate around conical intersection, because || ▽H || -> 0 at asymptote

    !dim x NStates x NStates 3-order tensor ▽H
    !Transform ▽H to (▽H)^2 representation
    !eigval harvests eigen values of (▽H)^2
    !eigvec harvests eigen vectors of (▽H)^2 in representation same to input ▽H
    subroutine NondegenerateRepresentation(dH,eigval,eigvec,dim,NStates)
        integer,intent(in)::dim,NStates
        real*8,dimension(dim,NStates,NStates),intent(inout)::dH
        real*8,dimension(NStates),intent(out)::eigval
        real*8,dimension(NStates,NStates),intent(out)::eigvec
        logical::degenerate
        integer::i,j,k,l
        eigvec=sy3matdotmul(dH,dH,dim,NStates)
        call My_dsyev('V',eigvec,eigval,NStates)
        dH=sy3UnitaryTransformation(dH,eigvec,dim,NStates)
        call CheckDegeneracy(degenerate,AlmostDegenerate,eigval,NStates)
        if(degenerate) write(*,*)'Warning: nondegenerate representation is also almost degenerate'
    end subroutine NondegenerateRepresentation
!------------------------- End --------------------------

end module Basic