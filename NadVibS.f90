!Generate NadVibS input files
module NadVibS
    use Basic
    use ElectronicStructure
    use Analyzation
    implicit none

!Example: type(HdExpansionCoefficient),allocatable,dimension(:,:)::HdEC
    !         HdEC(jstate,istate).Order(iorder).Array(i) is the i-th expansion coefficient
    !         in iorder-th order terms for Hd(jstate,istate)
    type HdExpansionCoefficient
        type(d2PArray),allocatable,dimension(:)::order
    end type HdExpansionCoefficient

contains
subroutine GenerateNadVibSInput()
    integer::ip,istate,jstate,iorder
    real*8,dimension(InternalDimension)::qPrecursor,qSuccessor,freqPrecursor,freqSuccessor
    real*8,dimension(InternalDimension,InternalDimension)::HPrecursor,HSuccessor
    real*8,dimension(InternalDimension,CartesianDimension)::BPrecursor,BSuccessor
    real*8,dimension(InternalDimension,InternalDimension,NState,NState)::Htemp
    type(Data),allocatable,dimension(:)::pointtemp
    !Definition of dshift and Tshift see Schuurman & Yarkony 2008 JCP 128 eq. (12)
    real*8,dimension(InternalDimension)::dshift
    real*8,dimension(InternalDimension,InternalDimension)::Tshift
    !Precursor
    call WilsonBMatrixAndInternalCoordinateq(BPrecursor,qPrecursor,reshape(MoleculeDetail.RefConfig,[CartesianDimension]),InternalDimension,CartesianDimension)
    call ReadElectronicStructureHessian(HPrecursor,InternalDimension)
    call WilsonGFMethod(freqPrecursor,HPrecursor,InternalDimension,BPrecursor,MoleculeDetail.mass,MoleculeDetail.NAtoms)
    if(minval(freqPrecursor)<-1d-14) write(*,*)'Warning: imaginary frequency found for precursor'
    !Successor
        !Allocate storage space
        allocate(pointtemp(NPoints))
        do ip=1,NPoints
            allocate(pointtemp(ip).geom(CartesianDimension))
            allocate(pointtemp(ip).energy(NState))
            allocate(pointtemp(ip).dH(CartesianDimension,NState,NState))
        end do
    call ReadElectronicStructureData(pointtemp,NPoints)!Read reference geometry
    call WilsonBMatrixAndInternalCoordinateq(BSuccessor,qSuccessor,pointtemp(IndexReference).geom,InternalDimension,CartesianDimension)
!This version directly use ground state Hessian at reference geometry
Htemp=AdiabaticddH(qSuccessor)
HSuccessor=Htemp(:,:,1,1)
!I will write a version shift the reference to ground state minimum of Hd someday
    call WilsonGFMethod(freqSuccessor,HSuccessor,InternalDimension,BSuccessor,MoleculeDetail.mass,MoleculeDetail.NAtoms)
    if(minval(freqSuccessor)<-1d-14) write(*,*)'Warning: imaginary frequency found for successor'
    dshift=matmul(transpose(HPrecursor),qSuccessor-qPrecursor)
    Tshift=matmul(transpose(HPrecursor),HSuccessor)
    open(unit=99,file='nadvibs.in',status='replace')
        write(99,'(A59)')'Angular frequency in atomic unit of each vibrational basis:'
        write(99,*)freqSuccessor
        do istate=1,NState
            do jstate=istate,NState
                do iorder=0,NOrder
                    write(99,'(A2,I2,I2,A15,I2)')'Hd',jstate,istate,'Expansion order',iorder
                    write(99,*)Hd_HdEC(jstate,istate).Order(iorder).Array
                end do
            end do
        end do
        write(99,'(A63)')'Angular frequency in atomic unit of each precursor normal mode:'
        write(99,*)freqPrecursor
        write(99,'(A13)')'Shift Vector:'
        write(99,*)dshift
        write(99,'(A22)')'Transformation Matrix:'
        write(99,*)Tshift
    close(99)
end subroutine GenerateNadVibSInput

end module NadVibs