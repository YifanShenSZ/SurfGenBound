!Diabatic Hamiltonian (Hd) and related quantities
module DiabaticHamiltonian
    use General
    use Mathematics
    use LinearAlgebra
    implicit none

!Derived type
    !Example: type(HdExpansionCoefficient),allocatable,dimension(:,:)::HdEC
    !         HdEC(jstate,istate).Order(iorder).Array(i) is the i-th expansion coefficient
    !         in iorder-th order terms for Hd(jstate,istate)
    type HdExpansionCoefficient
        type(d2PArray),allocatable,dimension(:)::order
    end type HdExpansionCoefficient

    !This defines the mapping rule between an actual basis and its serial number
    !Example: type(ExpansionBasisDefinition),allocatable,dimension(:)::EBNR
    !         EBNR(i).order is the polynomial order for i-th basis function
    !         This function is the product of q(EBNR(i).indice(1))*...*q(EBNR(i).indice(EBNR(i).order))
    !         Where q is the internal coordinate (for this program it's internal coordinate difference)
    type ExpansionBasisNumberingRule
        !In this program, expansion basis function is polynomial of internal coordinate
        integer::order!The order of this polynomial
        integer,allocatable,dimension(:)::indice!The product of which internal coordinates forms this polynomial
    end type ExpansionBasisNumberingRule

!Global variable
    !Number of Hd expansion basis functions, number of expansion coefficients
    integer::NExpansionBasis,NExpansionCoefficients

!DiabaticHamiltonian module only variable
    !Basic information of Hd:
    !    NStates: number of states Hd describes
    !    InternalDimension: dimension of internal space Hd takes into account
    !    NOrder: Hd expansion order
    !    HdEC & EBNR: see derived type section above
    integer::Hd_NStates,Hd_InternalDimension,Hd_NOrder
    type(HdExpansionCoefficient),allocatable,dimension(:,:)::Hd_HdEC!Short for Hd Expansion Coefficient, use only lower triangle
    type(ExpansionBasisNumberingRule),allocatable,dimension(:)::Hd_EBNR!short for Expansion Basis Numbering Rule

contains
!The initializer for DiabaticHamiltonian module
subroutine InitializeDiabaticHamiltonian(NewHd,NStates,NOrder)
    logical,intent(in),optional::NewHd
    integer,intent(in),optional::NStates,NOrder
    logical::degenerate
    integer::istate,jstate,iorder,i,n
    !Initialize Hd expansion coefficient (Hd_HdEC)
        if((.not.present(NewHd)).or.(.not.NewHd)) then
            fill0=NewHd
        else
            fill0=.false.
        end if
        select case(JobType)
            case('FitNewDiabaticHamiltonian')!To fit Hd from scratch, provide an initial guess
                !Allocate storage space
                    allocate(Hd_HdEC(Hd_NStates,Hd_NStates))
                    do istate=1,Hd_NStates
                        do jstate=istate,Hd_NStates
                            allocate(Hd_HdEC(jstate,istate).Order(0:Hd_NOrder))
                            do iorder=0,Hd_NOrder
                                allocate(Hd_HdEC(jstate,istate).Order(iorder).Array(int(iCombination(Hd_InternalDimension+iorder-1,iorder))))
                                Hd_HdEC(jstate,istate).Order(iorder).Array=0d0
                            end do
                        end do
                    end do
                call CheckDegeneracy(degenerate,AlmostDegenerate,ReferencePoint.energy,Hd_NStates)
                if(Degenerate) then
                    forall(istate=1:Nstates,jstate=1:Hd_NStates,istate>=jstate)
                        Hd_HdEC(istate,jstate).Order(0).Array(1)=ReferencePoint.H(istate,jstate)
                    end forall
                else
                    forall(istate=2:Hd_NStates)
                        Hd_HdEC(istate,istate).Order(0).Array(1)=ReferencePoint.energy(istate)-ReferencePoint.energy(1)
                    end forall
                end if
                forall(istate=1:Nstates,jstate=1:Hd_NStates,istate>=jstate)
                    Hd_HdEC(istate,jstate).Order(1).Array=ReferencePoint.dH(:,istate,jstate)
                end forall
            case('ContinueFitting')!Read old Hd expansion coefficients
                !Allocate storage space
                    allocate(Hd_HdEC(Hd_NStates,Hd_NStates))
                    do istate=1,Hd_NStates
                        do jstate=istate,Hd_NStates
                            allocate(Hd_HdEC(jstate,istate).Order(0:Hd_NOrder))
                            do iorder=0,Hd_NOrder
                                allocate(Hd_HdEC(jstate,istate).Order(iorder).Array(int(iCombination(Hd_InternalDimension+iorder-1,iorder))))
                                Hd_HdEC(jstate,istate).Order(iorder).Array=0d0
                            end do
                        end do
                    end do
                call ReadHdExpansionCoefficients()
                if(ReferenceChange) stop 'Program abort: reference point changed, not supported yet'
            case default!Use same NStates & NOrder of the fitted Hd, read old Hd expansion coefficients
                open(unit=99,file='HdExpansionCoefficient.out',status='old')
                    read(99,*)
                    read(99,*)Hd_NStates
                    read(99,*)
                    read(99,*)Hd_NOrder
                close(99)
                !Allocate storage space
                    allocate(Hd_HdEC(Hd_NStates,Hd_NStates))
                    do istate=1,Hd_NStates
                        do jstate=istate,Hd_NStates
                            allocate(Hd_HdEC(jstate,istate).Order(0:Hd_NOrder))
                            do iorder=0,Hd_NOrder
                                allocate(Hd_HdEC(jstate,istate).Order(iorder).Array(int(iCombination(Hd_InternalDimension+iorder-1,iorder))))
                                Hd_HdEC(jstate,istate).Order(iorder).Array=0d0
                            end do
                        end do
                    end do
                call ReadHdExpansionCoefficients()
        end select
    !Initialize NExpansionBasis, NExpansionCoefficients
        NExpansionBasis=0!Set counter to 0
        !Count how many basis functions for an Hd element
        do i=0,Hd_NOrder
            NExpansionBasis=NExpansionBasis+size(Hd_HdEC(1,1).Order(i).Array)
        end do
        !Times the number of independent Hd elements, we obtain the number of expansion coefficients
        NExpansionCoefficients=Hd_NStates*(Hd_NStates+1)/2*NExpansionBasis
    call InitializeExpansionBasisNumberingRule()
end subroutine InitializeDiabaticHamiltonian

!Load Hd expansion coefficient from HdExpansionCoefficient.out to global variable Hd_HdEC
subroutine ReadHdExpansionCoefficients()
    integer::NStates,NOrder!The old Hd is not necessarily fitted under same condition
    integer::istate,jstate,iorder
    open(unit=99,file='HdExpansionCoefficient.out',status='old')
        read(99,*)
        read(99,*)NStates
        read(99,*)
        read(99,*)NOrder
        read(99,*)
        do istate=1,NStates
            do jstate=istate,NStates
                read(99,*)
                do iorder=0,NOrder
                    read(99,*)
                    read(99,*)Hd_HdEC(jstate,istate).Order(iorder).Array
                end do
            end do
        end do
    close(99)
end subroutine ReadHdExpansionCoefficients

!Write the specified Hd expansion coefficient to file HdExpansionCoefficient.out 
subroutine WriteHdExpansionCoefficients(HdEC,NStates,NOrder)
    type(HdExpansionCoefficient),dimension(NStates,NStates),intent(in)::HdEC
    integer,intent(in)::NStates,NOrder
    integer::istate,jstate,iorder
    open(unit=99,file='HdExpansionCoefficient.out',status='replace')
        write(99,'(A28)')'Number of electronic states:'
        write(99,*)NStates
        write(99,'(A37)')'Diabatic hamiltonian expansion order:'
        write(99,*)NOrder
        write(99,'(A40)')'Expansion coefficients for each element:'
        do istate=1,NStates
            do jstate=istate,NStates
                write(99,'(A2,I2,I2)')'Hd',jstate,istate
                do iorder=0,NOrder
                    write(99,'(A15,I2)')'Expansion order',iorder
                    write(99,*)Hd_HdEC(jstate,istate).Order(iorder).Array
                end do
            end do
        end do
    close(99)
end subroutine WriteHdExpansionCoefficients

!-------------- Hd definition ---------------
    !This version generates Hd in NadVibS format: cast the elements of Hd into
    !simple polynomials of the internal coordinate difference from reference geometry
    !
    !Comment:
    !This treatment is good only for bounded system, where the molecule is semi-rigid

    subroutine InitializeExpansionBasisNumberingRule()!Hd_EBNR
        integer::i,iorder,n
        i=0!Count the number of the expansion basis functions
        do iorder=0,Hd_NOrder
            i=i+int(iCombination(Hd_InternalDimension+iorder-1,iorder))
        end do
        allocate(Hd_EBNR(i))
        n=1!The serial number of the expansion basis function
        do iorder=0,Hd_NOrder
            !The 1st one in each order
            Hd_EBNR(n).order=iorder
            allocate(Hd_EBNR(n).indice(iorder))
            Hd_EBNR(n).indice=1
            n=n+1
            !This is a pseudo counter, we add 1 to the 1st digit then carry to latter digits
            !but it should satisfy Hd_EBNR(n).indice(i)>=Hd_EBNR(n).indice(i+1)
            do istate=2,int(iCombination(Hd_InternalDimension+iorder-1,iorder))
                Hd_EBNR(n).order=iorder
                allocate(Hd_EBNR(n).indice(iorder))
                Hd_EBNR(n).indice=Hd_EBNR(n-1).indice
                !Add 1 to the 1st digit
                Hd_EBNR(n).indice(1)=Hd_EBNR(n).indice(1)+1
                !Carry to latter digits
                do i=1,iorder
                    if(Hd_EBNR(n).indice(i)>Hd_InternalDimension) then
                        Hd_EBNR(n).indice(i)=1
                        Hd_EBNR(n).indice(i+1)=Hd_EBNR(n).indice(i+1)+1
                    end if
                end do
                !Modify to satisfy Hd_EBNR(n).indice(i)>=Hd_EBNR(n).indice(i+1)
                do i=iorder-1,1,-1
                    if(Hd_EBNR(n).indice(i)<Hd_EBNR(n).indice(i+1)) then
                        Hd_EBNR(n).indice(i)=Hd_EBNR(n).indice(i+1)
                    end if
                end do
                n=n+1
            end do
        end do
    end subroutine InitializeExpansionBasisNumberingRule

    !The value of n-th expansion basis function at some coordinate q
    function ExpansionBasis(q,n)
        real*8::ExpansionBasis
        real*8,dimension(Hd_InternalDimension),intent(in)::q
        integer,intent(in)::n
        integer::i
        ExpansionBasis=1d0
        do i=1,Hd_EBNR(n).order
            ExpansionBasis=ExpansionBasis*q(Hd_EBNR(n).indice(i))
        end do
    end function ExpansionBasis

    !The value of ▽(n-th expansion basis function) at some coordinate q
    function ExpansionBasisGradient(q,n)
        real*8,dimension(Hd_InternalDimension)::ExpansionBasisGradient
        real*8,dimension(Hd_InternalDimension),intent(in)::q
        integer,intent(in)::n
        integer::m,i,OrderCount
        do m=1,Hd_InternalDimension
            OrderCount=0
            do i=1,Hd_EBNR(n).order
                if (Hd_EBNR(n).indice(i)==m) OrderCount=OrderCount+1
            end do
            if(OrderCount>0) then
                ExpansionBasisGradient(m)=dble(OrderCount)*q(m)**(OrderCount-1)
                do i=1,Hd_EBNR(n).order
                    if(Hd_EBNR(n).indice(i)/=m) ExpansionBasisGradient(m)=ExpansionBasisGradient(m)*q(Hd_EBNR(n).indice(i))
                end do
            else
                ExpansionBasisGradient(m)=0d0
            end if
        end do
    end function ExpansionBasisGradient

    !The value of ▽▽(n-th expansion basis function) at some coordinate q
    function ExpansionBasisHessian(q,n)
        real*8,dimension(Hd_InternalDimension,Hd_InternalDimension)::ExpansionBasisHessian
        real*8,dimension(Hd_InternalDimension),intent(in)::q
        integer,intent(in)::n
        integer::m1,m2,i,OrderCount1,OrderCount2
        do m1=1,Hd_InternalDimension
            OrderCount1=0!Diagonal
            do i=1,Hd_EBNR(n).order
                if (Hd_EBNR(n).indice(i)==m1) OrderCount1=OrderCount1+1
            end do
            if(OrderCount1>0) then
                ExpansionBasisHessian(m1,m1)=dble(OrderCount1*(OrderCount1-1))*q(m1)**(OrderCount1-2)
                do i=1,Hd_EBNR(n).order
                    if(Hd_EBNR(n).indice(i)/=m1) ExpansionBasisHessian(m1,m1)=ExpansionBasisHessian(m1,m1)*q(Hd_EBNR(n).indice(i))
                end do
            else
                ExpansionBasisHessian(m1,m1)=0d0
            end if
            do m2=m1+1,Hd_InternalDimension!Off-diagonal
                OrderCount1=0
                OrderCount2=0
                do i=1,Hd_EBNR(n).order
                    if(Hd_EBNR(n).indice(i)==m1) OrderCount1=OrderCount1+1
                    if(Hd_EBNR(n).indice(i)==m2) OrderCount2=OrderCount2+1
                end do
                if(OrderCount1>0.and.OrderCount2>0) then
                    ExpansionBasisHessian(m2,m1)=dble(OrderCount1*OrderCount2)*q(m1)**(OrderCount1-1)*q(m2)**(OrderCount2-1)
                    do i=1,Hd_EBNR(n).order
                        if(Hd_EBNR(n).indice(i)/=m1.and.Hd_EBNR(n).indice(i)/=m2) ExpansionBasisHessian(m2,m1)=ExpansionBasisHessian(m2,m1)*q(Hd_EBNR(n).indice(i))
                    end do
                else
                    ExpansionBasisHessian(m2,m1)=0d0
                end if
            end do
        end do
    end function ExpansionBasisHessian
!------------------- End --------------------

!------------ Diabatic quantity -------------
    function Hd(q)!Return the value of Hd in diabatic representation at some coordinate q
        real*8,dimension(Hd_NStates,Hd_NStates)::Hd
        real*8,dimension(Hd_InternalDimension),intent(in)::q
        integer::istate,jstate,iorder,i,n
        real*8,dimension(NExpansionBasis)::f
        do i=1,NExpansionBasis
            f(i)=ExpansionBasis(q,i)
        end do
        do istate=1,Hd_NStates
            do jstate=istate,Hd_NStates
                Hd(jstate,istate)=0d0
                n=1!The serial number of the expansion basis function
                do iorder=0,Hd_NOrder
                    do i=1,size(Hd_HdEC(jstate,istate).Order(iorder).Array)
                        Hd(jstate,istate)=Hd(jstate,istate)+Hd_HdEC(jstate,istate).Order(iorder).Array(i)*f(n)
                        n=n+1
                    end do
                end do
            end do
        end do
    end function Hd

    function dHd(q)!Return the value of ▽Hd in diabatic representation at some coordinate q
        real*8,dimension(Hd_InternalDimension,Hd_NStates,Hd_NStates)::dHd
        real*8,dimension(Hd_InternalDimension),intent(in)::q
        integer::istate,jstate,iorder,i,n
        real*8,dimension(Hd_InternalDimension,NExpansionBasis)::fd
        do i=1,NExpansionBasis
            fd(:,i)=ExpansionBasisGradient(q,i)
        end do
        do istate=1,Hd_NStates
            do jstate=istate,Hd_NStates
                dHd(:,jstate,istate)=0d0
                n=1!The serial number of the expansion basis function
                do iorder=0,Hd_NOrder
                    do i=1,size(Hd_HdEC(jstate,istate).Order(iorder).Array)
                        dHd(:,jstate,istate)=dHd(:,jstate,istate)&
                            +Hd_HdEC(jstate,istate).Order(iorder).Array(i)*fd(:,n)
                        n=n+1
                    end do
                end do
            end do
        end do
    end function dHd

    function ddHd(q)!Return the value of ▽▽Hd in diabatic representation at some coordinate q
        real*8,dimension(Hd_InternalDimension,Hd_InternalDimension,Hd_NStates,Hd_NStates)::ddHd
        real*8,dimension(Hd_InternalDimension),intent(in)::q
        integer::istate,jstate,iorder,i,n
        real*8,dimension(Hd_InternalDimension,Hd_InternalDimension,NExpansionBasis)::fdd
        do i=1,NExpansionBasis
            fdd(:,:,i)=ExpansionBasisHessian(q,i)
        end do
        do istate=1,Hd_NStates
            do jstate=istate,Hd_NStates
                ddHd(:,:,jstate,istate)=0d0
                n=1!The serial number of the expansion basis function
                do iorder=0,Hd_NOrder
                    do i=1,size(Hd_HdEC(jstate,istate).Order(iorder).Array)
                        ddHd(:,:,jstate,istate)=ddHd(:,:,jstate,istate)&
                            +Hd_HdEC(jstate,istate).Order(iorder).Array(i)*fdd(:,:,n)
                        n=n+1
                    end do
                end do
            end do
        end do
    end function ddHd

    !The value of Hd in diabatic representation and expansion basis functions at some coordinate q
    subroutine Hd_f(Hd,f,q)
        real*8,dimension(Hd_NStates,Hd_NStates),intent(out)::Hd
        real*8,dimension(NExpansionBasis),intent(out)::f
        real*8,dimension(Hd_InternalDimension),intent(in)::q
        integer::istate,jstate,iorder,i,n
        do i=1,NExpansionBasis
            f(i)=ExpansionBasis(q,i)
        end do
        do istate=1,Hd_NStates
            do jstate=istate,Hd_NStates
                Hd(jstate,istate)=0d0
                n=1!The serial number of the expansion basis function
                do iorder=0,Hd_NOrder
                    do i=1,size(Hd_HdEC(jstate,istate).Order(iorder).Array)
                        Hd(jstate,istate)=Hd(jstate,istate)&
                            +Hd_HdEC(jstate,istate).Order(iorder).Array(i)*f(n)
                        n=n+1
                    end do
                end do
            end do
        end do
    end subroutine Hd_f

    !The value of ▽Hd in diabatic representation and expansion basis function gradient at some coordinate q
    !fd(:,i) stores the gradient of i-th expansion basis function
    subroutine dHd_fd(dHd,fd,q)
        real*8,dimension(Hd_InternalDimension,Hd_NStates,Hd_NStates),intent(out)::dHd
        real*8,dimension(Hd_InternalDimension,NExpansionBasis),intent(out)::fd
        real*8,dimension(Hd_InternalDimension),intent(in)::q
        integer::istate,jstate,iorder,i,n
        do i=1,NExpansionBasis
            fd(:,i)=ExpansionBasisGradient(q,i)
        end do
        do istate=1,Hd_NStates
            do jstate=istate,Hd_NStates
                dHd(:,jstate,istate)=0d0
                n=1!The serial number of the expansion basis function
                do iorder=0,Hd_NOrder
                    do i=1,size(Hd_HdEC(jstate,istate).Order(iorder).Array)
                        dHd(:,jstate,istate)=dHd(:,jstate,istate)&
                            +Hd_HdEC(jstate,istate).Order(iorder).Array(i)*fd(:,n)
                        n=n+1
                    end do
                end do
            end do
        end do
    end subroutine dHd_fd
!------------------- End --------------------

!------------ Adiabatic quantity ------------
    !Compute adiabatic quantity from Hd at some coordinate q

    function AdiabaticEnergy(q)!Return adiabatic energy
        real*8,dimension(Hd_NStates)::AdiabaticEnergy
        real*8,dimension(Hd_InternalDimension),intent(in)::q
        real*8,dimension(Hd_NStates,Hd_NStates)::phi
        phi=Hd(q)
        call My_dsyev('N',phi,AdiabaticEnergy,Hd_NStates)
    end function AdiabaticEnergy

    function AdiabaticdH(q)!Return adiabatic gradient
        real*8,dimension(Hd_InternalDimension,Hd_NStates,Hd_NStates)::AdiabaticdH
        real*8,dimension(Hd_InternalDimension),intent(in)::q
        real*8,dimension(Hd_NStates)::energy
        real*8,dimension(Hd_NStates,Hd_NStates)::phi
        phi=Hd(q)
        call My_dsyev('V',phi,energy,Hd_NStates)
        AdiabaticdH=dHd(q)
        AdiabaticdH=sy3UnitaryTransformation(AdiabaticdH,phi,Hd_InternalDimension,Hd_NStates)
    end function AdiabaticdH

    function AdiabaticddH(q)!Return adiabatic Hessian
        real*8,dimension(Hd_InternalDimension,Hd_InternalDimension,Hd_NStates,Hd_NStates)::AdiabaticddH
        real*8,dimension(Hd_InternalDimension),intent(in)::q
        real*8,dimension(Hd_InternalDimension,Hd_NStates,Hd_NStates)::dH,M
        real*8,dimension(Hd_NStates)::energy
        real*8,dimension(Hd_NStates,Hd_NStates)::phi
        phi=Hd(q)
        call My_dsyev('V',phi,energy,Hd_NStates)
        dH=dHd(q)
        dH=sy3UnitaryTransformation(dH,phi,Hd_InternalDimension,Hd_NStates)
        M=deigvec_ByKnowneigval_dA(energy,dH,Hd_InternalDimension,Hd_NStates)
        AdiabaticddH=asy3matdirectmulsy3(M,dH,Hd_InternalDimension,Hd_InternalDimension,Hd_NStates)
        AdiabaticddH=sy4UnitaryTransformation(ddHd(q),phi,Hd_InternalDimension,Hd_InternalDimension,Hd_NStates)&
            -AdiabaticddH-transpose4(AdiabaticddH,Hd_InternalDimension,Hd_InternalDimension,Hd_NStates,Hd_NStates)
    end function AdiabaticddH

    !energy harvests adiabatic energy, dH harvests ▽H_a
    subroutine AdiabaticEnergy_dH(q,energy,dH)
        real*8,dimension(Hd_InternalDimension),intent(in)::q
        real*8,dimension(Hd_NStates),intent(out)::energy
        real*8,dimension(Hd_InternalDimension,Hd_NStates,Hd_NStates),intent(out)::dH
        real*8,dimension(Hd_NStates,Hd_NStates)::phi
        phi=Hd(q)
        call My_dsyev('V',phi,energy,Hd_NStates)
        dH=dHd(q)
        dH=sy3UnitaryTransformation(dH,phi,Hd_InternalDimension,Hd_NStates)
    end subroutine AdiabaticEnergy_dH

    !phi harvests adiabatic states in diabatic representation, f harvests expansion basis function values 
    subroutine AdiabaticEnergy_State_f(q,energy,phi,f)
        real*8,dimension(Hd_InternalDimension),intent(in)::q
        real*8,dimension(Hd_NStates),intent(out)::energy
        real*8,dimension(Hd_NStates,Hd_NStates),intent(out)::phi
        real*8,dimension(NExpansionBasis),intent(out)::f
        call Hd_f(phi,f,q)
        call My_dsyev('V',phi,energy,Hd_NStates)
    end subroutine AdiabaticEnergy_State_f

    !fd harvests expansion basis function gradient values, fd(:,i) = the gradient of i-th expansion basis function
    subroutine AdiabaticEnergy_dH_State_f_fd(q,energy,dH,phi,f,fd)
        real*8,dimension(Hd_InternalDimension),intent(in)::q
        real*8,dimension(Hd_NStates),intent(out)::energy
        real*8,dimension(Hd_InternalDimension,Hd_NStates,Hd_NStates),intent(out)::dH
        real*8,dimension(Hd_NStates,Hd_NStates),intent(out)::phi
        real*8,dimension(NExpansionBasis),intent(out)::f
        real*8,dimension(Hd_InternalDimension,NExpansionBasis),intent(out)::fd
        call Hd_f(phi,f,q)
        call My_dsyev('V',phi,energy,Hd_NStates)
        call dHd_fd(dH,fd,q)
        dH=sy3UnitaryTransformation(dH,phi,Hd_InternalDimension,Hd_NStates)
    end subroutine AdiabaticEnergy_dH_State_f_fd
!------------------- End --------------------

!---------- Nondegenerate quantity ----------
    !When Hamiltonian is almost degenerate, adiabatic basis experiences numerical difficulty
    !This causes no trouble to H, but operators other than H suffer a lot (e.g. ▽H)
    !Diabatz is always good, but it cannot be determined without knowledge of neighbourhood,
    !sometimes we ask for a nondegenerate representation based only on single point information
    !Choose (▽H)^2 representation: it is good enough for fitting a C1 molecule when not exploring too much seam space
    !    (▽H)^2 is Hermitian so there exists exactly one (▽H)^2 representation
    !    It reduces to s orthogonal to h for 2-fold degeneracy case, where s is the average force of 2 PESs
    !Note it is only appropriate around conical intersection, because || ▽H || -> 0 at asymptote
    !Compute nondegenerate representation quantity from Hd at some coordinate q

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

    !H harvests H_nd, dH harvests ▽H_nd
    subroutine NondegenerateH_dH(q,H,dH)
        real*8,dimension(Hd_InternalDimension),intent(in)::q
        real*8,dimension(Hd_NStates,Hd_NStates),intent(out)::H
        real*8,dimension(Hd_InternalDimension,Hd_NStates,Hd_NStates),intent(out)::dH
        real*8,dimension(Hd_NStates)::eigval
        real*8,dimension(Hd_NStates,Hd_NStates)::phi
        H=Hd(q)
        dH=dHd(q)
        call NondegenerateRepresentation(dH,eigval,phi,Hd_InternalDimension,Hd_NStates)
        H=matmul(transpose(phi),matmul(H,phi))
    end subroutine NondegenerateH_dH

    !phi harvests the nondegenerate states in diabatic representation
    !f harvests expansion basis function values
    !fd harvests expansion basis function gradient values, fd(:,i) = the gradient of i-th expansion basis function
    subroutine NondegenerateH_dH_State_f_fd(q,H,dH,phi,f,fd)
        real*8,dimension(Hd_InternalDimension),intent(in)::q
        real*8,dimension(Hd_NStates,Hd_NStates),intent(out)::H
        real*8,dimension(Hd_InternalDimension,Hd_NStates,Hd_NStates),intent(out)::dH
        real*8,dimension(Hd_NStates,Hd_NStates),intent(out)::phi
        real*8,dimension(NExpansionBasis),intent(out)::f
        real*8,dimension(Hd_InternalDimension,NExpansionBasis),intent(out)::fd
        real*8,dimension(Hd_NStates)::eigval
        call Hd_f(H,f,q)
        call dHd_fd(dH,fd,q)
        call NondegenerateRepresentation(dH,eigval,phi,Hd_InternalDimension,Hd_NStates)
        H=matmul(transpose(phi),matmul(H,phi))
    end subroutine NondegenerateH_dH_State_f_fd

    !eigval harvests the eigenvalues corresponding to phi
    !dHd harvests ▽H_d
    subroutine NondegenerateH_dH_eigval_State_dHd_f_fd(q,H,dH,eigval,phi,dHd,f,fd)
        real*8,dimension(Hd_InternalDimension),intent(in)::q
        real*8,dimension(Hd_NStates,Hd_NStates),intent(out)::H
        real*8,dimension(Hd_InternalDimension,Hd_NStates,Hd_NStates),intent(out)::dH
        real*8,dimension(Hd_NStates),intent(out)::eigval
        real*8,dimension(Hd_NStates,Hd_NStates),intent(out)::phi
        real*8,dimension(Hd_InternalDimension,Hd_NStates,Hd_NStates),intent(out)::dHd
        real*8,dimension(NExpansionBasis),intent(out)::f
        real*8,dimension(Hd_InternalDimension,NExpansionBasis),intent(out)::fd
        call Hd_f(H,f,q)
        call dHd_fd(dHd,fd,q)
        dH=dHd
        call NondegenerateRepresentation(dH,eigval,phi,Hd_InternalDimension,Hd_NStates)
        H=matmul(transpose(phi),matmul(H,phi))
    end subroutine NondegenerateH_dH_eigval_State_dHd_f_fd
!------------------- End --------------------

end module DiabaticHamiltonian