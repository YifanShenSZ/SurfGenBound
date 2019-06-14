module NadVibSInterface
    use Basic
    use ElectronicStructure
    implicit none

!Derived type
    !Example: type(NVS_HdExpansionCoefficient),allocatable,dimension(:,:)::NVS_HdEC
    !         NVS_HdEC(jstate,istate).Order(iorder).Array(i) is the i-th expansion coefficient
    !         in iorder-th order terms for Hd(jstate,istate)
    type NVS_HdExpansionCoefficient!Store Hd expansion coefficient in NadVibS format
        type(d2PArray),allocatable,dimension(:)::Order
    end type NVS_HdExpansionCoefficient

!NadVibSInterface module only variable
    integer::NVS_NOrder,NVS_NHdExpansionBasis
    integer,allocatable,dimension(:)::NVS_NumberOfEachOrderTerms
    type(ExpansionBasisNumberingRule),allocatable,dimension(:)::NVS_EBNR
    type(NVS_HdExpansionCoefficient),allocatable,dimension(:,:)::NVS_HdEC

contains
subroutine GenerateNadVibSInput()
    character*2::chartemp
    integer::i,j,k
    real*8::dbletemp
    real*8,dimension(InternalDimension)::qPrecursor,qSuccessor,freqPrecursor,freqSuccessor
    real*8,dimension(CartesianDimension)::rSuccesor
    real*8,dimension(InternalDimension,InternalDimension)::HPrecursor,HSuccessor,modePrecursor,modeSuccessor
    real*8,dimension(InternalDimension,CartesianDimension)::BPrecursor,BSuccessor
    real*8,dimension(InternalDimension,InternalDimension,NState,NState)::Htemp
    !Definition of dshift and Tshift see Schuurman & Yarkony 2008 JCP 128 eq. (12)
    real*8,dimension(InternalDimension)::dshift
    real*8,dimension(InternalDimension,InternalDimension)::Tshift
    call InitializeNadVibSInterface()
    !Precursor
    call WilsonBMatrixAndInternalCoordinateq(BPrecursor,qPrecursor,reshape(MoleculeDetail.RefConfig,[CartesianDimension]),InternalDimension,CartesianDimension)
    call ReadElectronicStructureHessian(HPrecursor,InternalDimension)
    call WilsonGFMethod(freqPrecursor,modePrecursor,HPrecursor,InternalDimension,BPrecursor,MoleculeDetail.mass,MoleculeDetail.NAtoms)
    if(minval(freqPrecursor)<0d0) stop 'Program abort: imaginary frequency found for precursor'
    !Successor
!This version chooses the neutral ground state minimum as origin, corresponding normal mode as basis
    open(unit=99,file='MinimumCartesianGeometry.xyz',status='old')
    	read(99,*)
    	read(99,*)
        do i=1,MoleculeDetail.NAtoms
            read(99,'(A2,3F20.15)')chartemp,rSuccesor(3*i-2:3*i)
        end do
        rSuccesor=rSuccesor*AInAU!Convert to atomic unit
    close(99)
    call WilsonBMatrixAndInternalCoordinateq(BSuccessor,qSuccessor,rSuccesor,InternalDimension,CartesianDimension)
    qSuccessor=qSuccessor-ReferencePoint.geom
    Htemp=AdiabaticddH(qSuccessor)
    HSuccessor=Htemp(:,:,1,1)
!End of the specific treatment
    call WilsonGFMethod(freqSuccessor,modeSuccessor,HSuccessor,InternalDimension,BSuccessor,MoleculeDetail.mass,MoleculeDetail.NAtoms)
    if(minval(freqSuccessor)<0d0) stop 'Program abort: imaginary frequency found for successor'
    write(*,'(1x,A53)')'Suggestion on number of basis by distance estimation:'
    dshift=dAbs(matmul(modeSuccessor,qSuccessor-qPrecursor))
    do i=1,InternalDimension
        dshift(i)=dshift(i)/(2d0*freqSuccessor(i))
        do j=3,7,2
            dbletemp=dFactorial2(j-2)
            if(dbletemp**(1d0/dble(j-1))>dshift(i)) exit
        end do
        if(j==3.and.freqSuccessor(i)/cm_1InAU>3000d0) j=1!Short distance X-H stretching usually can be frozen
        write(*,'(5x,A4,I3,A14,I2)')'Mode',i,', Basis number',j
    end do
    call OriginShift(qSuccessor)!Shift origin to ground state minimum
    call HdEc_Hd2NVS(modeSuccessor)!Reformat Hd expansion coefficient into NadVibS format
    !Definition of dshift and Tshift see Schuurman & Yarkony 2008 JCP 128 eq. (12)
    dshift=matmul(modePrecursor,qSuccessor-qPrecursor)
    call My_dgetri(modeSuccessor,InternalDimension)
    Tshift=matmul(modePrecursor,modeSuccessor)
    open(unit=99,file='nadvibs.in',status='replace')
        write(99,'(A54)')'Angular frequency of each vibrational basis: (In a.u.)'
        write(99,*)freqSuccessor
        do i=1,NState
            do j=i,NState
                do k=0,NVS_NOrder
                    write(99,'(A2,I2,I2,A15,I2)')'Hd',j,i,'Expansion order',k
                    write(99,*)NVS_HdEC(j,i).Order(k).Array
                end do
            end do
        end do
        write(99,'(A58)')'Angular frequency of each precursor normal mode: (In a.u.)'
        write(99,*)freqPrecursor
        write(99,'(A13)')'Shift Vector:'
        write(99,*)dshift
        write(99,'(A22)')'Transformation Matrix:'
        write(99,*)Tshift
    close(99)
end subroutine GenerateNadVibSInput

subroutine InitializeNadVibSInterface()
    integer::i,j,iorder,n
    NVS_NOrder=Hd_EBNR(1).order!Determine the highest order used
    NVS_NHdExpansionBasis=0!Count the number of the expansion basis functions
    allocate(NVS_NumberOfEachOrderTerms(0:NVS_NOrder))
    do iorder=0,NVS_NOrder
        NVS_NumberOfEachOrderTerms(iorder)=iCombination(InternalDimension+iorder-1,iorder)
        NVS_NHdExpansionBasis=NVS_NHdExpansionBasis+NVS_NumberOfEachOrderTerms(iorder)
    end do
    allocate(NVS_EBNR(NVS_NHdExpansionBasis))!Generate expansion basis numbering mapping (NVS_EBNR)
    n=1!The serial number of the expansion basis function
    do iorder=0,NVS_NOrder
        NVS_EBNR(n).order=iorder!The 1st one in each order
        allocate(NVS_EBNR(n).indice(iorder))
        NVS_EBNR(n).indice=1
        n=n+1
        !This is a pseudo counter, we add 1 to the 1st digit then carry to latter digits
        !but it should satisfy NVS_EBNR(n).indice(i)>=NVS_EBNR(n).indice(i+1)
        do j=2,NVS_NumberOfEachOrderTerms(iorder)
            NVS_EBNR(n).order=iorder
            allocate(NVS_EBNR(n).indice(iorder))
            NVS_EBNR(n).indice=NVS_EBNR(n-1).indice
            !Add 1 to the 1st digit
            NVS_EBNR(n).indice(1)=NVS_EBNR(n).indice(1)+1
            !Carry to latter digits
            do i=1,iorder
                if(NVS_EBNR(n).indice(i)>InternalDimension) then
                    NVS_EBNR(n).indice(i)=1
                    NVS_EBNR(n).indice(i+1)=NVS_EBNR(n).indice(i+1)+1
                end if
            end do
            !Modify to satisfy NVS_EBNR(n).indice(i)>=NVS_EBNR(n).indice(i+1)
            do i=iorder-1,1,-1
                if(NVS_EBNR(n).indice(i)<NVS_EBNR(n).indice(i+1)) then
                    NVS_EBNR(n).indice(i)=NVS_EBNR(n).indice(i+1)
                end if
            end do
            n=n+1
        end do
    end do
    allocate(NVS_HdEC(NState,NState))!Allocate storage space of NVS_HdEC
    do j=1,NState
        do i=j,NState
            allocate(NVS_HdEC(i,j).Order(0:NVS_NOrder))
            do iorder=0,NVS_NOrder
                allocate(NVS_HdEC(i,j).Order(iorder).Array(NVS_NumberOfEachOrderTerms(iorder)))
                NVS_HdEC(i,j).Order(iorder).Array=0d0
            end do
        end do
    end do
end subroutine InitializeNadVibSInterface

subroutine HdEc_Hd2NVS(mode)
    real*8,dimension(InternalDimension,InternalDimension),intent(in)::mode
    !do iorder=0,NVS_NOrder!Fill in the order by order form
    !    indice=1
    !    i=WhichExpansionBasis(iorder,indice(1:iorder))
    !    forall(istate=1:NState,jstate=1:NState,istate>=jstate)
    !        HdEC(istate,jstate).Order(iorder).Array(1)=Hd_HdEC(istate,jstate).Array(i)
    !    end forall
    !    do j=2,size(HdEC(1,1).Order(iorder).Array)
    !        indice(1)=indice(1)+1
    !        do i=1,iorder
    !            if(indice(i)>InternalDimension) then
    !                indice(i)=1
    !                indice(i+1)=indice(i+1)+1
    !            end if
    !        end do
    !        do i=iorder-1,1,-1
    !            if(indice(i)<indice(i+1)) indice(i)=indice(i+1)
    !        end do
    !        i=WhichExpansionBasis(iorder,indice(1:iorder))
    !        forall(istate=1:NState,jstate=1:NState,istate>=jstate)
    !            HdEC(istate,jstate).Order(iorder).Array(j)=Hd_HdEC(istate,jstate).Array(i)
    !        end forall
    !    end do
    !end do
end subroutine HdEc_Hd2NVS

integer function NVS_WhichExpansionBasis(order,indice)
    integer,intent(in)::order
    integer,dimension(order),intent(inout)::indice
    integer,dimension(order)::temp
    integer::low,up
    low=0
    do up=0,order-1
        low=low+NVS_NumberOfEachOrderTerms(up)
    end do
    up=low+NVS_NumberOfEachOrderTerms(order)
    low=low+1
    call iQuickSort(indice,1,order,temp,order)
    call bisect(low,up)
    contains
    recursive subroutine bisect(low,up)
        integer,intent(in)::low,up
        integer::i,bisection
        if(low==up) then
            NVS_WhichExpansionBasis=low
        else if(up-low==1) then
            do i=order,1,-1
                if(indice(i)/=NVS_EBNR(low).indice(i)) exit
            end do
            if(i<1) then
                NVS_WhichExpansionBasis=low
            else
                NVS_WhichExpansionBasis=up
            end if
        else
            bisection=(low+up)/2
            do i=order,1,-1
                if(indice(i)/=NVS_EBNR(bisection).indice(i)) exit
            end do
            if(i<1) then
                NVS_WhichExpansionBasis=bisection
            else
                if(indice(i)>NVS_EBNR(bisection).indice(i)) then
                    call bisect(bisection,up)
                else
                    call bisect(low,bisection)
                end if
            end if
        end if
    end subroutine bisect
end function NVS_WhichExpansionBasis

end module NadVibSInterface