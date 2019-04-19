!This version generates Hd in NadVibs format: cast the elements of Hd into
!simple polynomials of the internal coordinate difference from the reference geometry
!
!Comment:
!This treatment is good only for bounded system, where the molecule is semi-rigid
module DiabaticHamiltonian
    use Basic
    implicit none

!Derived type
    !Example: type(HdExpansionCoefficient),allocatable,dimension(:,:)::HdEC
    !         HdEC(jstate,istate).Order(iorder).Array(i) is the i-th expansion coefficient
    !         in iorder-th order terms for HdEC(jstate,istate)
    type HdExpansionCoefficient
        type(d2PArray),allocatable,dimension(:)::order
    end type HdExpansionCoefficient

    !To form linear equations, we have to number the expansion basis functions
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

!DiabaticHamiltonian global variable
    !Number of Hd expansion basis functions, number of expansion coefficients
    integer::NExpansionBasis,NExpansionCoefficients
    type(HdExpansionCoefficient),allocatable,dimension(:,:)::HdEC!Short for Hd Expansion Coefficient, use only lower triangle
    type(ExpansionBasisNumberingRule),allocatable,dimension(:)::EBNR!short for Expansion Basis Numbering Rule

contains
!The initializer for DiabaticHamiltonian module
subroutine InitializeDiabaticHamiltonian()
    logical::degenerate
    integer::istate,jstate,iorder,i,n
    !Initialize Hd expansion coefficient (HdEC)
        select case(JobType)
            case('FitNewDiabaticHamiltonian')!To fit Hd from scratch, provide an initial guess
                !Allocate storage space
                    allocate(HdEC(NStates,NStates))
                    do istate=1,NStates
                        do jstate=istate,NStates
                            allocate(HdEC(jstate,istate).Order(0:NOrder))
                            do iorder=0,NOrder
                                allocate(HdEC(jstate,istate).Order(iorder).Array(int(iCombination(InternalDimension+iorder-1,iorder))))
                                HdEC(jstate,istate).Order(iorder).Array=0d0
                            end do
                        end do
                    end do
                call CheckDegeneracy(degenerate,AlmostDegenerate,ReferencePoint.energy,NStates)
                if(Degenerate) then
                    forall(istate=1:Nstates,jstate=1:NStates,istate>=jstate)
                        HdEC(istate,jstate).Order(0).Array(1)=ReferencePoint.H(istate,jstate)
                    end forall
                else
                    forall(istate=2:NStates)
                        HdEC(istate,istate).Order(0).Array(1)=ReferencePoint.energy(istate)-ReferencePoint.energy(1)
                    end forall
                end if
                forall(istate=1:Nstates,jstate=1:NStates,istate>=jstate)
                    HdEC(istate,jstate).Order(1).Array=ReferencePoint.dH(:,istate,jstate)
                end forall
            case('ContinueFitting')!Read old Hd expansion coefficients
                !Allocate storage space
                    allocate(HdEC(NStates,NStates))
                    do istate=1,NStates
                        do jstate=istate,NStates
                            allocate(HdEC(jstate,istate).Order(0:NOrder))
                            do iorder=0,NOrder
                                allocate(HdEC(jstate,istate).Order(iorder).Array(int(iCombination(InternalDimension+iorder-1,iorder))))
                                HdEC(jstate,istate).Order(iorder).Array=0d0
                            end do
                        end do
                    end do
                call ReadHdExpansionCoefficients()
                if(ReferenceChange) stop 'Program abort: reference point changed, not supported yet'
            case('NadVibS')!Use same NStates & NOrder of the fitted Hd, read old Hd expansion coefficients
                open(unit=99,file='HdExpansionCoefficient.out',status='old')
                    read(99,*)
                    read(99,*)NStates
                    read(99,*)
                    read(99,*)NOrder
                close(99)
                !Allocate storage space
                    allocate(HdEC(NStates,NStates))
                    do istate=1,NStates
                        do jstate=istate,NStates
                            allocate(HdEC(jstate,istate).Order(0:NOrder))
                            do iorder=0,NOrder
                                allocate(HdEC(jstate,istate).Order(iorder).Array(int(iCombination(InternalDimension+iorder-1,iorder))))
                                HdEC(jstate,istate).Order(iorder).Array=0d0
                            end do
                        end do
                    end do
                call ReadHdExpansionCoefficients()
            case default
        end select
    !Initialize NExpansionBasis, NExpansionCoefficients
        NExpansionBasis=0!Set counter to 0
        !Count how many basis functions for an Hd element
        do i=0,NOrder
            NExpansionBasis=NExpansionBasis+size(HdEC(1,1).Order(i).Array)
        end do
        !Times the number of independent Hd elements, we obtain the number of expansion coefficients
        NExpansionCoefficients=NStates*(NStates+1)/2*NExpansionBasis
    !Initialize expansion basis numbering mapping (EBNR)
        i=0!Count the number of the expansion basis functions
        do iorder=0,NOrder
            i=i+int(iCombination(InternalDimension+iorder-1,iorder))
        end do
        allocate(EBNR(i))
        n=1!The serial number of the expansion basis function
        do iorder=0,NOrder
            !The 1st one in each order
            EBNR(n).order=iorder
            allocate(EBNR(n).indice(iorder))
            EBNR(n).indice=1
            n=n+1
            !This is a pseudo counter, we add 1 to the 1st digit then carry to latter digits
            !but it should satisfy EBNR(n).indice(i)>=EBNR(n).indice(i+1)
            do istate=2,int(iCombination(InternalDimension+iorder-1,iorder))
                EBNR(n).order=iorder
                allocate(EBNR(n).indice(iorder))
                EBNR(n).indice=EBNR(n-1).indice
                !Add 1 to the 1st digit
                EBNR(n).indice(1)=EBNR(n).indice(1)+1
                !Carry to latter digits
                do i=1,iorder
                    if(EBNR(n).indice(i)>InternalDimension) then
                        EBNR(n).indice(i)=1
                        EBNR(n).indice(i+1)=EBNR(n).indice(i+1)+1
                    end if
                end do
                !Modify to satisfy EBNR(n).indice(i)>=EBNR(n).indice(i+1)
                do i=iorder-1,1,-1
                    if(EBNR(n).indice(i)<EBNR(n).indice(i+1)) then
                        EBNR(n).indice(i)=EBNR(n).indice(i+1)
                    end if
                end do
                n=n+1
            end do
        end do
end subroutine InitializeDiabaticHamiltonian

!Load Hd expansion coefficient from HdExpansionCoefficient.out to global variable HdEC
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
                    read(99,*)HdEC(jstate,istate).Order(iorder).Array
                end do
            end do
        end do
    close(99)
end subroutine ReadHdExpansionCoefficients

!Write the specified Hd expansion coefficient to file HdExpansionCoefficient.out 
subroutine WriteHdExpansionCoefficients(HdEC,NStates,NOrders)
    type(HdExpansionCoefficient),dimension(NStates,NStates),intent(in)::HdEC
    integer,intent(in)::NStates,NOrders
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
                    write(99,*)HdEC(jstate,istate).Order(iorder).Array
                end do
            end do
        end do
    close(99)
end subroutine WriteHdExpansionCoefficients

!The value of n-th expansion basis function at some coordinate q
function ExpansionBasis(q,n)
    real*8::ExpansionBasis
    real*8,dimension(InternalDimension),intent(in)::q
    integer,intent(in)::n
    integer::i
    ExpansionBasis=1d0
    do i=1,EBNR(n).order
        ExpansionBasis=ExpansionBasis*q(EBNR(n).indice(i))
    end do
end function ExpansionBasis

!The value of ▽(n-th expansion basis function) at some coordinate q
function ExpansionBasisGradient(q,n)
    real*8,dimension(InternalDimension)::ExpansionBasisGradient
    real*8,dimension(InternalDimension),intent(in)::q
    integer,intent(in)::n
    integer::m,i,OrderCount
    do m=1,InternalDimension
        OrderCount=0
        do i=1,EBNR(n).order
            if (EBNR(n).indice(i)==m) OrderCount=OrderCount+1
        end do
        if(OrderCount>0) then
            ExpansionBasisGradient(m)=dble(OrderCount)*q(m)**(OrderCount-1)
            do i=1,EBNR(n).order
                if(EBNR(n).indice(i)/=m) ExpansionBasisGradient(m)=ExpansionBasisGradient(m)*q(EBNR(n).indice(i))
            end do
        else
            ExpansionBasisGradient(m)=0d0
        end if
    end do
end function ExpansionBasisGradient

!The value of ▽▽(n-th expansion basis function) at some coordinate q
function ExpansionBasisHessian(q,n)
    real*8,dimension(InternalDimension,InternalDimension)::ExpansionBasisHessian
    real*8,dimension(InternalDimension),intent(in)::q
    integer,intent(in)::n
    integer::m1,m2,i,OrderCount1,OrderCount2
    do m1=1,InternalDimension
        OrderCount1=0!Diagonal
        do i=1,EBNR(n).order
            if (EBNR(n).indice(i)==m1) OrderCount1=OrderCount1+1
        end do
        if(OrderCount1>0) then
            ExpansionBasisHessian(m1,m1)=dble(OrderCount1*(OrderCount1-1))*q(m1)**(OrderCount1-2)
            do i=1,EBNR(n).order
                if(EBNR(n).indice(i)/=m1) ExpansionBasisHessian(m1,m1)=ExpansionBasisHessian(m1,m1)*q(EBNR(n).indice(i))
            end do
        else
            ExpansionBasisHessian(m1,m1)=0d0
        end if
        do m2=m1+1,InternalDimension!Off-diagonal
            OrderCount1=0
            OrderCount2=0
            do i=1,EBNR(n).order
                if(EBNR(n).indice(i)==m1) OrderCount1=OrderCount1+1
                if(EBNR(n).indice(i)==m2) OrderCount2=OrderCount2+1
            end do
            if(OrderCount1>0.and.OrderCount2>0) then
                ExpansionBasisHessian(m2,m1)=dble(OrderCount1*OrderCount2)*q(m1)**(OrderCount1-1)*q(m2)**(OrderCount2-1)
                do i=1,EBNR(n).order
                    if(EBNR(n).indice(i)/=m1.and.EBNR(n).indice(i)/=m2) ExpansionBasisHessian(m2,m1)=ExpansionBasisHessian(m2,m1)*q(EBNR(n).indice(i))
                end do
            else
                ExpansionBasisHessian(m2,m1)=0d0
            end if
        end do
    end do
end function ExpansionBasisHessian

!------------ Diabatic quantity -------------
    !The value of Hd in diabatic representation at some coordinate q
    function Hd(q)
        real*8,dimension(NStates,NStates)::Hd
        real*8,dimension(InternalDimension),intent(in)::q
        integer::istate,jstate,iorder,i,n
        real*8,dimension(NExpansionBasis)::f
        do i=1,NExpansionBasis
            f(i)=ExpansionBasis(q,i)
        end do
        do istate=1,NStates
            do jstate=istate,NStates
                Hd(jstate,istate)=0d0
                n=1!The serial number of the expansion basis function
                do iorder=0,NOrder
                    do i=1,size(HdEC(jstate,istate).Order(iorder).Array)
                        Hd(jstate,istate)=Hd(jstate,istate)+HdEC(jstate,istate).Order(iorder).Array(i)*f(n)
                        n=n+1
                    end do
                end do
            end do
        end do
    end function Hd

    !The value of Hd in diabatic representation and expansion basis functions at some coordinate q
    subroutine Hd_f(Hd,f,q)
        real*8,dimension(NStates,NStates),intent(out)::Hd
        real*8,dimension(NExpansionBasis),intent(out)::f
        real*8,dimension(InternalDimension),intent(in)::q
        integer::istate,jstate,iorder,i,n
        do i=1,NExpansionBasis
            f(i)=ExpansionBasis(q,i)
        end do
        do istate=1,NStates
            do jstate=istate,NStates
                Hd(jstate,istate)=0d0
                n=1!The serial number of the expansion basis function
                do iorder=0,NOrder
                    do i=1,size(HdEC(jstate,istate).Order(iorder).Array)
                        Hd(jstate,istate)=Hd(jstate,istate)&
                            +HdEC(jstate,istate).Order(iorder).Array(i)*f(n)
                        n=n+1
                    end do
                end do
            end do
        end do
    end subroutine Hd_f

    !The value of ▽Hd in diabatic representation at some coordinate q
    function dHd(q)
        real*8,dimension(InternalDimension,NStates,NStates)::dHd
        real*8,dimension(InternalDimension),intent(in)::q
        integer::istate,jstate,iorder,i,n
        real*8,dimension(InternalDimension,NExpansionBasis)::fd
        do i=1,NExpansionBasis
            fd(:,i)=ExpansionBasisGradient(q,i)
        end do
        do istate=1,NStates
            do jstate=istate,NStates
                dHd(:,jstate,istate)=0d0
                n=1!The serial number of the expansion basis function
                do iorder=0,NOrder
                    do i=1,size(HdEC(jstate,istate).Order(iorder).Array)
                        dHd(:,jstate,istate)=dHd(:,jstate,istate)&
                            +HdEC(jstate,istate).Order(iorder).Array(i)*fd(:,n)
                        n=n+1
                    end do
                end do
            end do
        end do
    end function dHd

    !The value of ▽Hd in diabatic representation and expansion basis function gradient at some coordinate q
    !fd(:,i) stores the gradient of i-th expansion basis function
    subroutine dHd_fd(dHd,fd,q)
        real*8,dimension(InternalDimension,NStates,NStates),intent(out)::dHd
        real*8,dimension(InternalDimension,NExpansionBasis),intent(out)::fd
        real*8,dimension(InternalDimension),intent(in)::q
        integer::istate,jstate,iorder,i,n
        do i=1,NExpansionBasis
            fd(:,i)=ExpansionBasisGradient(q,i)
        end do
        do istate=1,NStates
            do jstate=istate,NStates
                dHd(:,jstate,istate)=0d0
                n=1!The serial number of the expansion basis function
                do iorder=0,NOrder
                    do i=1,size(HdEC(jstate,istate).Order(iorder).Array)
                        dHd(:,jstate,istate)=dHd(:,jstate,istate)&
                            +HdEC(jstate,istate).Order(iorder).Array(i)*fd(:,n)
                        n=n+1
                    end do
                end do
            end do
        end do
    end subroutine dHd_fd

    !The value of ▽▽Hd in diabatic representation at some coordinate q
    function ddHd(q)
        real*8,dimension(InternalDimension,InternalDimension,NStates,NStates)::ddHd
        real*8,dimension(InternalDimension),intent(in)::q
        integer::istate,jstate,iorder,i,n
        real*8,dimension(InternalDimension,InternalDimension,NExpansionBasis)::fdd
        do i=1,NExpansionBasis
            fdd(:,:,i)=ExpansionBasisHessian(q,i)
        end do
        do istate=1,NStates
            do jstate=istate,NStates
                ddHd(:,:,jstate,istate)=0d0
                n=1!The serial number of the expansion basis function
                do iorder=0,NOrder
                    do i=1,size(HdEC(jstate,istate).Order(iorder).Array)
                        ddHd(:,:,jstate,istate)=ddHd(:,:,jstate,istate)&
                            +HdEC(jstate,istate).Order(iorder).Array(i)*fdd(:,:,n)
                        n=n+1
                    end do
                end do
            end do
        end do
    end function ddHd
!------------------- End --------------------

!------------ Adiabatic quantity ------------
    !Compute adiabatic quantity from Hd at some coordinate q

    !Return adiabatic energy
    function AdiabaticEnergy(q)
        real*8,dimension(NStates)::AdiabaticEnergy
        real*8,dimension(InternalDimension),intent(in)::q
        real*8,dimension(NStates,NStates)::phi
        phi=Hd(q)
        call My_dsyev('N',phi,AdiabaticEnergy,NStates)
    end function AdiabaticEnergy

    !Return adiabatic gradient
    function AdiabaticdH(q)
        real*8,dimension(InternalDimension,NStates,NStates)::AdiabaticdH
        real*8,dimension(InternalDimension),intent(in)::q
        real*8,dimension(NStates)::energy
        real*8,dimension(NStates,NStates)::phi
        phi=Hd(q)
        call My_dsyev('V',phi,energy,NStates)
        AdiabaticdH=dHd(q)
        AdiabaticdH=sy3UnitaryTransformation(AdiabaticdH,phi,InternalDimension,NStates)
    end function AdiabaticdH

    !Return adiabatic Hessian
    function AdiabaticddH(q)
        real*8,dimension(InternalDimension,InternalDimension,NStates,NStates)::AdiabaticddH
        real*8,dimension(InternalDimension),intent(in)::q
        real*8,dimension(InternalDimension,NStates,NStates)::dH,M
        real*8,dimension(NStates)::energy
        real*8,dimension(NStates,NStates)::phi
        phi=Hd(q)
        call My_dsyev('V',phi,energy,NStates)
        dH=dHd(q)
        dH=sy3UnitaryTransformation(dH,phi,InternalDimension,NStates)
        M=deigvec_ByKnowneigval_dA(energy,dH,InternalDimension,NStates)
        AdiabaticddH=asy3matdirectmulsy3(M,dH,InternalDimension,InternalDimension,NStates)
        AdiabaticddH=sy4UnitaryTransformation(ddHd(q),phi,InternalDimension,InternalDimension,NStates)&
            -AdiabaticddH-transpose4(AdiabaticddH,InternalDimension,InternalDimension,NStates,NStates)
    end function AdiabaticddH

    !energy harvests adiabatic energy, dH harvests ▽H_a
    subroutine AdiabaticEnergy_dH(q,energy,dH)
        real*8,dimension(InternalDimension),intent(in)::q
        real*8,dimension(NStates),intent(out)::energy
        real*8,dimension(InternalDimension,NStates,NStates),intent(out)::dH
        real*8,dimension(NStates,NStates)::phi
        phi=Hd(q)
        call My_dsyev('V',phi,energy,NStates)
        dH=dHd(q)
        dH=sy3UnitaryTransformation(dH,phi,InternalDimension,NStates)
    end subroutine AdiabaticEnergy_dH

    !phi harvests adiabatic states in diabatic representation, f harvests expansion basis function values 
    subroutine AdiabaticEnergy_State_f(q,energy,phi,f)
        real*8,dimension(InternalDimension),intent(in)::q
        real*8,dimension(NStates),intent(out)::energy
        real*8,dimension(NStates,NStates),intent(out)::phi
        real*8,dimension(NExpansionBasis),intent(out)::f
        call Hd_f(phi,f,q)
        call My_dsyev('V',phi,energy,NStates)
    end subroutine AdiabaticEnergy_State_f

    !fd harvests expansion basis function gradient values, fd(:,i) = the gradient of i-th expansion basis function
    subroutine AdiabaticEnergy_dH_State_f_fd(q,energy,dH,phi,f,fd)
        real*8,dimension(InternalDimension),intent(in)::q
        real*8,dimension(NStates),intent(out)::energy
        real*8,dimension(InternalDimension,NStates,NStates),intent(out)::dH
        real*8,dimension(NStates,NStates),intent(out)::phi
        real*8,dimension(NExpansionBasis),intent(out)::f
        real*8,dimension(InternalDimension,NExpansionBasis),intent(out)::fd
        call Hd_f(phi,f,q)
        call My_dsyev('V',phi,energy,NStates)
        call dHd_fd(dH,fd,q)
        dH=sy3UnitaryTransformation(dH,phi,InternalDimension,NStates)
    end subroutine AdiabaticEnergy_dH_State_f_fd
!------------------- End --------------------

!---------- Nondegenerate quantity ----------
    !Compute nondegenerate representation quantity from Hd at some coordinate q
    !For details of the nondegenerate operator O adopted here, see module Basic

    !H harvests H_nd, dH harvests ▽H_nd
    subroutine NondegenerateH_dH(q,H,dH)
        real*8,dimension(InternalDimension),intent(in)::q
        real*8,dimension(NStates,NStates),intent(out)::H
        real*8,dimension(InternalDimension,NStates,NStates),intent(out)::dH
        real*8,dimension(NStates)::eigval
        real*8,dimension(NStates,NStates)::phi
        H=Hd(q)
        dH=dHd(q)
        call NondegenerateRepresentation(dH,eigval,phi,InternalDimension,NStates)
        H=matmul(transpose(phi),matmul(H,phi))
    end subroutine NondegenerateH_dH

    !phi harvests the nondegenerate states in diabatic representation
    !f harvests expansion basis function values
    !fd harvests expansion basis function gradient values, fd(:,i) = the gradient of i-th expansion basis function
    subroutine NondegenerateH_dH_State_f_fd(q,H,dH,phi,f,fd)
        real*8,dimension(InternalDimension),intent(in)::q
        real*8,dimension(NStates,NStates),intent(out)::H
        real*8,dimension(InternalDimension,NStates,NStates),intent(out)::dH
        real*8,dimension(NStates,NStates),intent(out)::phi
        real*8,dimension(NExpansionBasis),intent(out)::f
        real*8,dimension(InternalDimension,NExpansionBasis),intent(out)::fd
        real*8,dimension(NStates)::eigval
        call Hd_f(H,f,q)
        call dHd_fd(dH,fd,q)
        call NondegenerateRepresentation(dH,eigval,phi,InternalDimension,NStates)
        H=matmul(transpose(phi),matmul(H,phi))
    end subroutine NondegenerateH_dH_State_f_fd

    !eigval harvests the eigenvalues corresponding to phi
    !dHd harvests ▽H_d
    subroutine NondegenerateH_dH_eigval_State_dHd_f_fd(q,H,dH,eigval,phi,dHd,f,fd)
        real*8,dimension(InternalDimension),intent(in)::q
        real*8,dimension(NStates,NStates),intent(out)::H
        real*8,dimension(InternalDimension,NStates,NStates),intent(out)::dH
        real*8,dimension(NStates),intent(out)::eigval
        real*8,dimension(NStates,NStates),intent(out)::phi
        real*8,dimension(InternalDimension,NStates,NStates),intent(out)::dHd
        real*8,dimension(NExpansionBasis),intent(out)::f
        real*8,dimension(InternalDimension,NExpansionBasis),intent(out)::fd
        call Hd_f(H,f,q)
        call dHd_fd(dHd,fd,q)
        dH=dHd
        call NondegenerateRepresentation(dH,eigval,phi,InternalDimension,NStates)
        H=matmul(transpose(phi),matmul(H,phi))
    end subroutine NondegenerateH_dH_eigval_State_dHd_f_fd
!------------------- End --------------------

end module DiabaticHamiltonian