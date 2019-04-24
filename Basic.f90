!Library interface, global parameter & derived type & data storage, basic routine
module Basic
    use General
    use Mathematics
    use LinearAlgebra
    use GeometryTransformation
    use NonlinearOptimization
    use Nonadiabatic
    use DiabaticHamiltonian
    implicit none

!Parameter
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

	type MoleculeDetails!Store the details of an example molecule
	    integer::NAtoms
        character*2,allocatable,dimension(:)::ElementSymbol
        real*8,allocatable,dimension(:,:)::RefConfig!Short for REFerence CONFIGuration
        real*8,allocatable,dimension(:)::mass
    end type MoleculeDetails
    
!Input variable
    !See input example for the meaning of each variable
    !Main input
        character*32::JobType
        integer::NState,NOrder
        logical::SameTrainingSet
        integer::NPoints,IndexReference
        character*128::ArtifactGeometryDataFile,ArtifactEnergyDataFile
        integer::NArtifactPoints
    type(MoleculeDetails)::MoleculeDetail!Molecule detail input

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
subroutine InitializeBasic()!Initialize Basic module and lower level libraries
    call InitializePhaseFixing(NState)
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
    real*8,dimension(NState)::eigval
    real*8,dimension(NState,NState)::eigvec
    !Count how many points are almost degenerate
    NDegpoints=0
    do ip=1,NPoints
        call CheckDegeneracy(degenerate,AlmostDegenerate,point(ip).energy,NState)
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
        call CheckDegeneracy(degenerate,AlmostDegenerate,point(ip).energy,NState)
        if(degenerate) then
            IndicesDeg(index)=ip
            degpoint(index).weight=point(ip).weight
            allocate(degpoint(index).geom(dim))
            degpoint(index).geom=point(ip).geom
            allocate(degpoint(index).energy(dim))
            degpoint(index).energy=point(ip).energy
            allocate(degpoint(index).dH(dim,NState,NState))
            degpoint(index).dH=point(ip).dH
            call NondegenerateRepresentation(degpoint(index).dH,eigval,eigvec,dim,NState,DegenerateThreshold=AlmostDegenerate)
            allocate(degpoint(index).H(NState,NState))
            degpoint(index).H=transpose(eigvec)
            forall(istate=1:NState)
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
                do istate=1,NState
                    do jstate=istate,NState
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
                do istate=1,NState
                    do jstate=istate,NState
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
                allocate(point(ip).energy(NState))
                read(99,*)point(ip).energy
                allocate(point(ip).H(NState,NState))
                read(99,*)point(ip).H
                allocate(point(ip).dH(dim,NState,NState))
                do istate=1,NState
                    do jstate=istate,NState
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
                do istate=1,NState
                    do jstate=istate,NState
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

end module Basic