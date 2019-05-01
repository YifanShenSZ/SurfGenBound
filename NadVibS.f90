!Generate NadVibS input files
module NadVibS
    use Basic
    use ElectronicStructure
    use Analyzation
    implicit none

!Derived type
    !Example: type(NadVibS_HdEC),allocatable,dimension(:,:)::HdEC
    !         HdEC(jstate,istate).Order(iorder).Array(i) is the i-th expansion coefficient
    !         in iorder-th order terms for Hd(jstate,istate)
    type NadVibS_HdEC!Store Hd expansion coefficient in NadVibS format
        type(d2PArray),allocatable,dimension(:)::order
    end type NadVibS_HdEC

contains
subroutine GenerateNadVibSInput()
    integer::i,j,istate,jstate,iorder,NOrder
    integer,allocatable,dimension(:)::indice
    real*8,dimension(InternalDimension)::qPrecursor,qSuccessor,freqPrecursor,freqSuccessor
    real*8,dimension(CartesianDimension)::rSuccesor
    real*8,dimension(InternalDimension,InternalDimension)::HPrecursor,HSuccessor
    real*8,dimension(InternalDimension,CartesianDimension)::BPrecursor,BSuccessor
    real*8,dimension(InternalDimension,InternalDimension,NState,NState)::Htemp
    type(Data),allocatable,dimension(:)::pointtemp
    type(NadVibS_HdEC),dimension(NState,NState)::HdEC
    !Definition of dshift and Tshift see Schuurman & Yarkony 2008 JCP 128 eq. (12)
    real*8,dimension(InternalDimension)::dshift
    real*8,dimension(InternalDimension,InternalDimension)::Tshift
    !Precursor
    call WilsonBMatrixAndInternalCoordinateq(BPrecursor,qPrecursor,reshape(MoleculeDetail.RefConfig,[CartesianDimension]),InternalDimension,CartesianDimension)
    call ReadElectronicStructureHessian(HPrecursor,InternalDimension)
    call WilsonGFMethod(freqPrecursor,HPrecursor,InternalDimension,BPrecursor,MoleculeDetail.mass,MoleculeDetail.NAtoms)
    if(minval(freqPrecursor)<-1d-14) write(*,*)'Warning: imaginary frequency found for precursor'
    !Successor
    Analyzation_state=1!Search for ground state minimum
    qSuccessor=ReferencePoint.geom
    call BFGS(AdiabaticEnergyInterface,AdiabaticGradientInterface,qSuccessor,InternalDimension,&
        fdd=AdiabaticHessianInterface,f_fd=AdiabaticEnergy_GradientInterface)
    allocate(pointtemp(NPoints))!Read reference Cartesian geometry to provide int2cart initial guess
    do i=1,NPoints
        allocate(pointtemp(i).geom(CartesianDimension))
        allocate(pointtemp(i).energy(NState))
        allocate(pointtemp(i).dH(CartesianDimension,NState,NState))
    end do
    call ReadElectronicStructureData(pointtemp,NPoints)
    rSuccesor=CartesianCoordinater(qSuccessor,CartesianDimension,InternalDimension,&
        mass=MoleculeDetail.mass,r0=pointtemp(IndexReference).geom)
    call WilsonBMatrixAndInternalCoordinateq(BSuccessor,qSuccessor,rSuccesor,InternalDimension,CartesianDimension)
    Htemp=AdiabaticddH(qSuccessor)
    HSuccessor=Htemp(:,:,1,1)
    call WilsonGFMethod(freqSuccessor,HSuccessor,InternalDimension,BSuccessor,MoleculeDetail.mass,MoleculeDetail.NAtoms)
    if(minval(freqSuccessor)<-1d-14) write(*,*)'Warning: imaginary frequency found for successor'
    !Reformat Hd expansion coefficient into NadVibS format
    call OriginShift(qSuccessor-ReferencePoint.geom)!Shift origin to ground state minimum
    NOrder=0!Determine the highest order used
    do i=1,NHdExpansionBasis
        if(Hd_EBNR(i).order>NOrder) NOrder=Hd_EBNR(i).order
    end do
    allocate(indice(NOrder))
    do jstate=1,NState!Allocate storage space of HdEC
        do istate=jstate,NState
            allocate(HdEC(istate,jstate).Order(0:NOrder))
            do iorder=0,NOrder
                allocate(HdEC(istate,jstate).Order(iorder).Array(iCombination(InternalDimension+iorder-1,iorder)))
                HdEC(istate,jstate).Order(iorder).Array=0d0
            end do
        end do
    end do
    do iorder=0,NOrder
        indice=1
        i=WhichExpansionBasis(iorder,indice(1:iorder))
        forall(istate=1:NState,jstate=1:NState,istate>=jstate)
            HdEC(istate,jstate).Order(iorder).Array(1)=Hd_HdEC(istate,jstate).Array(i)
        end forall
        do j=2,size(HdEC(1,1).Order(iorder).Array)
            indice(1)=indice(1)+1
            do i=1,iorder
                if(indice(i)>InternalDimension) then
                    indice(i)=1
                    indice(i+1)=indice(i+1)+1
                end if
            end do
            do i=iorder-1,1,-1
                if(indice(i)<indice(i+1)) indice(i)=indice(i+1)
            end do
            i=WhichExpansionBasis(iorder,indice(1:iorder))
            forall(istate=1:NState,jstate=1:NState,istate>=jstate)
                HdEC(istate,jstate).Order(iorder).Array(j)=Hd_HdEC(istate,jstate).Array(i)
            end forall
        end do
    end do
    !Definition of dshift and Tshift see Schuurman & Yarkony 2008 JCP 128 eq. (12)
    dshift=matmul(transpose(HPrecursor),qSuccessor-qPrecursor)
    Tshift=matmul(transpose(HPrecursor),HSuccessor)
    open(unit=99,file='nadvibs.in',status='replace')
        write(99,'(A59)')'Angular frequency in atomic unit of each vibrational basis:'
        write(99,*)freqSuccessor
        do istate=1,NState
            do jstate=istate,NState
                do iorder=0,NOrder
                    write(99,'(A2,I2,I2,A15,I2)')'Hd',jstate,istate,'Expansion order',iorder
                    write(99,*)HdEC(jstate,istate).Order(iorder).Array
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