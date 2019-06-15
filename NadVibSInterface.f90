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

    type NVS_ExpansionBasisNumberingRule
        type(i2PArray),allocatable,dimension(:)::Number
    end type NVS_ExpansionBasisNumberingRule

!NadVibSInterface module only variable
    integer::NVS_NOrder
    integer,allocatable,dimension(:)::NVS_NumberOfEachOrderTerms
    type(NVS_ExpansionBasisNumberingRule),allocatable,dimension(:)::NVS_EBNR
    type(NVS_HdExpansionCoefficient),allocatable,dimension(:,:)::NVS_HdEC

contains
subroutine GenerateNadVibSInput()
    !Everything required for nadvibs.in
    real*8,dimension(InternalDimension)::qPrecursor,qSuccessor,freqPrecursor,freqSuccessor
    real*8,dimension(CartesianDimension)::rSuccesor
    real*8,dimension(InternalDimension,InternalDimension)::modePrecursor,LPrecursor,HPrecursor,modeSuccessor,LSuccessor,HSuccessor
    real*8,dimension(InternalDimension,CartesianDimension)::BPrecursor,BSuccessor
    real*8,dimension(InternalDimension)::dshift
    real*8,dimension(InternalDimension,InternalDimension)::Tshift
    !Work space
    character*2::chartemp
    integer::i,j,k
    real*8::dbletemp1,dbletemp2
    real*8,dimension(NState)::energy
    real*8,dimension(InternalDimension)::qtemp
    real*8,dimension(InternalDimension,InternalDimension,NState,NState)::Htemp
    call InitializeNadVibSInterface()
    !Precursor
    call WilsonBMatrixAndInternalCoordinateq(BPrecursor,qPrecursor,reshape(MoleculeDetail.RefConfig,[CartesianDimension]),InternalDimension,CartesianDimension)
    call ReadElectronicStructureHessian(HPrecursor,InternalDimension)
    call WilsonGFMethod(freqPrecursor,modePrecursor,LPrecursor,HPrecursor,InternalDimension,BPrecursor,MoleculeDetail.mass,MoleculeDetail.NAtoms)
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
    qtemp=qSuccessor-ReferencePoint.geom
    Htemp=AdiabaticddH(qtemp)
    HSuccessor=Htemp(:,:,1,1)
!End of the specific treatment
    call WilsonGFMethod(freqSuccessor,modeSuccessor,LSuccessor,HSuccessor,InternalDimension,BSuccessor,MoleculeDetail.mass,MoleculeDetail.NAtoms)
    if(minval(freqSuccessor)<0d0) stop 'Program abort: imaginary frequency found for successor'
    write(*,'(1x,A64)')'Suggestion on number of basis by distance and energy estimation:'
    !largest standard deviation > precursor-successor distance
    !2nd highest energy < precursor-successor ground state energy difference < highest
    dshift=dAbs(matmul(modeSuccessor,qPrecursor-qSuccessor))
    energy=AdiabaticEnergy(qPrecursor-ReferencePoint.geom)-AdiabaticEnergy(qSuccessor-ReferencePoint.geom)
    dbletemp1=1d0
    dbletemp2=1d0
    do i=1,InternalDimension
        dshift(i)=dshift(i)*dSqrt(freqSuccessor(i))
        do j=1,9!Consider (j-1)-th excited state standard deviation
            if(dSqrt(dFactorial2(2*j-1)/2d0**j)>dshift(i)) exit
        end do
        k=ceiling(energy(1)/freqSuccessor(i)-0.5d0)+1
        write(*,'(5x,A4,I3,A19,I2,A3,I2)')'Mode',i,', Basis number from',j,' to',k
        dbletemp1=dbletemp1*dble(j)
        dbletemp2=dbletemp2*dble(max(j,k))
    end do
    write(*,*)'The total number of smallest basis is',dbletemp1
    if(dbletemp1>2d0**31d0-1d0) write(*,*)'Warning: this is',dbletemp1/(2d0**31d0-1d0),'times larger than 2^31 - 1'
    write(*,*)'The total number of  largest basis is',dbletemp2
    if(dbletemp2>2d0**31d0-1d0) write(*,*)'Warning: this is',dbletemp2/(2d0**31d0-1d0),'times larger than 2^31 - 1'
ENERGY=AdiabaticEnergy(qPrecursor-ReferencePoint.geom)
WRITE(*,*)ENERGY
    call OriginShift(qSuccessor-ReferencePoint.geom)!Shift origin to ground state minimum
ENERGY=AdiabaticEnergy(qPrecursor-qSuccessor)
WRITE(*,*)ENERGY
    call HdEC_Hd2NVS(LSuccessor)!Reformat Hd expansion coefficient into NadVibS format
ENERGY=NVS_AdiabaticEnergy(matmul(modeSuccessor,qPrecursor-qSuccessor))
WRITE(*,*)ENERGY
    !Definition of dshift and Tshift see Schuurman & Yarkony 2008 JCP 128 eq. (12)
    dshift=matmul(modePrecursor,qSuccessor-qPrecursor)
    Tshift=matmul(modePrecursor,LSuccessor)
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

FUNCTION NVS_ADIABATICENERGY(Q)
    real*8,dimension(InternalDimension),INTENT(IN)::q
    real*8,dimension(NState)::NVS_ADIABATICENERGY
    integer::i,j,k,iorder
    real*8::dbletemp
    real*8,dimension(NState,NSTATE)::H
    forall(i=1:NState,j=1:NState,i>=j)
        H(i,j)=NVS_HdEC(i,j).Order(0).Array(1)
    end forall
    do iorder=1,NVS_NOrder
        do k=1,NVS_NumberOfEachOrderTerms(iorder)
            dbletemp=NVS_ExpansionBasis(q,iorder,k)
            forall(i=1:NState,j=1:NState,i>=j)
                H(i,j)=H(i,j)+NVS_HdEC(i,j).Order(iorder).Array(k)*dbletemp
            end forall
        end do
    end do
    call My_dsyev('N',H,NVS_ADIABATICENERGY,NState)
    contains
    real*8 function NVS_ExpansionBasis(q,order,n)
        real*8,dimension(InternalDimension),intent(in)::q
        integer,intent(in)::order,n
        integer::i
        NVS_ExpansionBasis=1d0
        do i=1,order
            NVS_ExpansionBasis=NVS_ExpansionBasis*q(NVS_EBNR(order).Number(n).Array(i))
        end do
    end function NVS_ExpansionBasis
END FUNCTION NVS_ADIABATICENERGY

subroutine InitializeNadVibSInterface()
    integer::i,j,iorder
    NVS_NOrder=Hd_EBNR(1).order!Determine the highest order used
    allocate(NVS_NumberOfEachOrderTerms(0:NVS_NOrder))
    allocate(NVS_EBNR(0:NVS_NOrder))
    do iorder=0,NVS_NOrder
        NVS_NumberOfEachOrderTerms(iorder)=iCombination(InternalDimension+iorder-1,iorder)
        allocate(NVS_EBNR(iorder).Number(NVS_NumberOfEachOrderTerms(iorder)))
        !Generate expansion basis numbering mapping for iorder-th order (NVS_EBNR(iorder))
        allocate(NVS_EBNR(iorder).Number(1).Array(iorder))
        NVS_EBNR(iorder).Number(1).Array=1
        !This is a pseudo InternalDimension+1 counter, former digit >= latter digit
        do j=2,NVS_NumberOfEachOrderTerms(iorder)
            allocate(NVS_EBNR(iorder).Number(j).Array(iorder))
            NVS_EBNR(iorder).Number(j).Array=NVS_EBNR(iorder).Number(j-1).Array
            NVS_EBNR(iorder).Number(j).Array(1)=NVS_EBNR(iorder).Number(j).Array(1)+1!Add 1 to the 1st digit
            do i=1,iorder-1!Carry to latter digits
                if(NVS_EBNR(iorder).Number(j).Array(i)>InternalDimension) then
                    NVS_EBNR(iorder).Number(j).Array(i)=1
                    NVS_EBNR(iorder).Number(j).Array(i+1)=NVS_EBNR(iorder).Number(j).Array(i+1)+1
                end if
            end do
            do i=iorder-1,1,-1!Modify to satisfy former digit >= latter digit
                if(NVS_EBNR(iorder).Number(j).Array(i)<NVS_EBNR(iorder).Number(j).Array(i+1)) &
                    NVS_EBNR(iorder).Number(j).Array(i)=NVS_EBNR(iorder).Number(j).Array(i+1)
            end do
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

subroutine HdEC_Hd2NVS(L)!Transform from internal coordinate Hd_HdEC to normal mode NVS_HdEC
    !This is done by:
    !    1, select an Hd_HdEC term
    !    2, under rotation, the terms making up the multiplication become the linear combination of every mode,
    !       i.e. q = L . Q, so we go through all combination and add the contribution to NVS_HdEC
    !    3, go to 1 until all Hd_HdEC are done
    real*8,dimension(InternalDimension,InternalDimension),intent(in)::L
    integer::order,location,n,i,j
    integer,dimension(NVS_NOrder)::indice,Hdindice
    real*8::coeff
    do n=1,NHdExpansionBasis!Main loop
        order=Hd_EBNR(n).order
        if(order==0) then!Const term will not change under any transformation
            forall(i=1:NState,j=1:NState,i>=j)
                NVS_HdEC(i,j).Order(0).Array(1)=NVS_HdEC(i,j).Order(0).Array(1)+Hd_HdEC(i,j).Array(n)
            end forall
        else!Go through all combination in a InternalDimension+1 counter manner
            Hdindice(1:order)=Hd_EBNR(n).indice
            indice(1:order)=1!indice(i)=j means using j-th mode at i-th position in multiplication
            do while(indice(order)<=InternalDimension)!Done when the counter overflows
                coeff=1d0
                do i=1,order
                    coeff=coeff*L(Hdindice(i),indice(i))
                end do
                if(coeff/=0d0) then
                    location=NVS_WhichExpansionBasis(order,indice(1:order))
                    forall(i=1:NState,j=1:NState,i>=j)
                        NVS_HdEC(i,j).Order(order).Array(location)=NVS_HdEC(i,j).Order(order).Array(location)&
                            +coeff*Hd_HdEC(i,j).Array(n)
                    end forall
                end if
                indice(1)=indice(1)+1!Add 1 to the InternalDimension+1 counter
                do i=1,order-1
                    if(indice(i)>InternalDimension) then!Carry
                        indice(i)=1
                        indice(i+1)=indice(i+1)+1
                    end if
                end do
            end do
        end if
    end do
end subroutine HdEC_Hd2NVS

integer function NVS_WhichExpansionBasis(order,indiceinput)
    integer,intent(in)::order
    integer,dimension(order),intent(in)::indiceinput
    integer,dimension(order)::indice,temp
    indice=-indiceinput
    call iQuickSort(indice,1,order,temp,order)!I only coded ascending sort
    indice=-indice
    call bisect(1,NVS_NumberOfEachOrderTerms(order))
    contains
    recursive subroutine bisect(low,up)
        integer,intent(in)::low,up
        integer::i,bisection
        if(low==up) then
            NVS_WhichExpansionBasis=low
        else if(up-low==1) then
            do i=order,1,-1
                if(indice(i)/=NVS_EBNR(order).Number(low).Array(i)) exit
            end do
            if(i<1) then
                NVS_WhichExpansionBasis=low
            else
                NVS_WhichExpansionBasis=up
            end if
        else
            bisection=(low+up)/2
            do i=order,1,-1
                if(indice(i)/=NVS_EBNR(order).Number(bisection).Array(i)) exit
            end do
            if(i<1) then
                NVS_WhichExpansionBasis=bisection
            else
                if(indice(i)>NVS_EBNR(order).Number(bisection).Array(i)) then
                    call bisect(bisection,up)
                else
                    call bisect(low,bisection)
                end if
            end if
        end if
    end subroutine bisect
end function NVS_WhichExpansionBasis

end module NadVibSInterface