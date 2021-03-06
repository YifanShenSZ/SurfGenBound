!Diabatic Hamiltonian (Hd) and related quantities
!EBNR and routines in Hd definition section are the only things to modify for basis change
module DiabaticHamiltonian
    use General; use Mathematics; use LinearAlgebra
	use Chemistry
    implicit none

!Derived type
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
    integer::NHdExpansionBasis,NHdExpansionCoefficients

!DiabaticHamiltonian module only variable
    !Basic information of Hd:
    !    NState: number of states Hd describes
    !    intdim: dimension of internal space Hd takes into account
    !    EBNR: see derived type section above
    !    HdEC: HdEC(i,j).Array(k) is the expansion coefficient for k-th basis of Hd(i,j)
    integer::Hd_NState,Hd_intdim
    type(ExpansionBasisNumberingRule),allocatable,dimension(:)::Hd_EBNR!short for Expansion Basis Numbering Rule
    type(d2PArray),allocatable,dimension(:,:)::Hd_HdEC!Short for Hd Expansion Coefficient, use only lower triangle

contains
!Initialize DiabaticHamiltonian module
!Required: NState & intdim
!Optional: NewHd: (default = false) if true, will return a blank Hd expansion coefficient
subroutine InitializeDiabaticHamiltonian(NState,intdim,NewHd)
    integer,intent(in)::NState,intdim
    logical,intent(in),optional::NewHd
    integer::istate,jstate,iorder,i,n
    Hd_NState=NState; Hd_intdim=intdim
    call InitializeExpansionBasisNumberingRule()!NHdExpansionBasis, EBNR
    NHdExpansionCoefficients=Hd_NState*(Hd_NState+1)/2*NHdExpansionBasis!Multiply the number of independent Hd elements
	allocate(Hd_HdEC(Hd_NState,Hd_NState))!HdEC
    do istate=1,Hd_NState
        do jstate=istate,Hd_NState
            allocate(Hd_HdEC(jstate,istate).Array(NHdExpansionBasis))
            Hd_HdEC(jstate,istate).Array=0d0
        end do
	end do
    if((.not.present(NewHd)).or.(.not.NewHd)) call ReadHdExpansionCoefficients(Hd_HdEC)
end subroutine InitializeDiabaticHamiltonian

!-------------- Hd definition ---------------
    !This version generates Hd in NadVibS format: cast the elements of Hd into
    !simple polynomials of the internal coordinate difference from reference geometry
    !
    !Comment:
    !This treatment is good only for bounded system, where the molecule is semi-rigid

    subroutine InitializeExpansionBasisNumberingRule()!Hd_EBNR
        integer::i,j
        NHdExpansionBasis=0
        open(unit=99,file='basis.in',status='old')
            do!Count how many expansion basis functions there are
                read(99,*,iostat=i); if(i/=0) exit
                NHdExpansionBasis=NHdExpansionBasis+1
            end do
            rewind 99
            allocate(Hd_EBNR(NHdExpansionBasis))!Read the definition of expansion basis functions
            do i=1,NHdExpansionBasis
                read(99,'(I5)',advance='no')Hd_EBNR(i).order
                allocate(Hd_EBNR(i).indice(Hd_EBNR(i).order))
                if(Hd_EBNR(i).order>0) then
                    do j=1,Hd_EBNR(i).order-1
                        read(99,'(I5)',advance='no')Hd_EBNR(i).indice(j)
                    end do
                    read(99,'(I5)')Hd_EBNR(i).indice(Hd_EBNR(i).order)
                else
                    read(99,*)
                end if
            end do
        close(99)
    end subroutine InitializeExpansionBasisNumberingRule

    integer function WhichExpansionBasis(order,indice)!Return the location of the specified basis in Hd_EBNR, 0 if not found
        integer,intent(in)::order
        integer,dimension(order),intent(in)::indice
        call bisect(1,NHdExpansionBasis)
        contains
        recursive subroutine bisect(low,up)
            integer,intent(in)::low,up
            integer::bisection,i
            if(up-low==1) then
                if(order==Hd_EBNR(low).order) then
                    do i=order,1,-1
                        if(indice(i)/=Hd_EBNR(low).indice(i)) exit
                    end do
                    if(i<1) then
                        WhichExpansionBasis=low
                        return
                    end if
                end if
                if(order==Hd_EBNR(up).order) then
                    do i=order,1,-1
                        if(indice(i)/=Hd_EBNR(up).indice(i)) exit
                    end do
                    if(i<1) then
                        WhichExpansionBasis=up
                        return
                    end if
                end if
                WhichExpansionBasis=0
            else
                bisection=(low+up)/2
                if(order>Hd_EBNR(bisection).order) then
                    call bisect(low,bisection)
                else if(order<Hd_EBNR(bisection).order) then
                    call bisect(bisection,up)
                else
                    do i=order,1,-1
                        if(indice(i)/=Hd_EBNR(bisection).indice(i)) exit
                    end do
                    if(i<1) then
                        WhichExpansionBasis=bisection
                    else
                        if(indice(i)>Hd_EBNR(bisection).indice(i)) then
                            call bisect(bisection,up)
                        else
                            call bisect(low,bisection)
                        end if
                    end if
                end if
            end if
        end subroutine bisect
    end function WhichExpansionBasis

    subroutine OriginShift(shift)!Transform HdEC according to origin shift from q0 to q1: shift = q1 - q0
        !This is done by:
        !    1, select an Hd expansion basis (Hd_HdEC)
        !    2, under translation, the terms making up the multiplication become var + const, so
        !       we go through all combination of var & const and add the contribution to HdECtemp
        !    3, go to 1 until all Hd_HdEC are done
        !    4, copy HdECtemp to Hd_HdEC
        real*8,dimension(Hd_intdim),intent(in)::shift
        integer::location,nvar,nconst,ivar,iconst,n,i,j
        integer,dimension(Hd_EBNR(1).order)::usevar,indicevar,indiceconst
        real*8::coeff
        type(d2PArray),dimension(Hd_NState,Hd_Nstate)::HdECtemp
        do j=1,Hd_NState!Allocate work space
            do i=j,Hd_NState
                allocate(HdECtemp(i,j).Array(NHdExpansionBasis))
                HdECtemp(i,j).Array=0d0
            end do
        end do
        do n=1,NHdExpansionBasis!Main loop
            if(Hd_EBNR(n).order==0) then!Const term will not change under any transformation
                forall(i=1:Hd_NState,j=1:Hd_NState,i>=j)
                    HdECtemp(i,j).Array(n)=HdECtemp(i,j).Array(n)+Hd_HdEC(i,j).Array(n)
                end forall
            else!Go through all combination in a binary counter manner
                usevar(1:Hd_EBNR(n).order)=0!usevar(i)=1 means using const at i-th position in multiplication
                do while(usevar(Hd_EBNR(n).order)<2)!Done when the counter overflows
                    nvar=sum(usevar(1:Hd_EBNR(n).order))
                    nconst=Hd_EBNR(n).order-nvar
                    ivar=1
                    iconst=1
                    do i=1,Hd_EBNR(n).order
                        if(usevar(i)==1) then
                            indicevar(ivar)=Hd_EBNR(n).indice(i)
                            ivar=ivar+1
                        else
                            indiceconst(iconst)=Hd_EBNR(n).indice(i)
                            iconst=iconst+1
                        end if
                    end do
                    location=WhichExpansionBasis(nvar,indicevar(1:nvar))
                    if(location==0) stop 'Program abort: basis space is not closed under origin shift'
                    coeff=1d0
                    do i=1,nconst
                        coeff=coeff*shift(indiceconst(i))
                    end do
                    forall(i=1:Hd_NState,j=1:Hd_NState,i>=j)
                        HdECtemp(i,j).Array(location)=HdECtemp(i,j).Array(location)+coeff*Hd_HdEC(i,j).Array(n)
                    end forall
                    usevar(1)=usevar(1)+1!Add 1 to the binary counter
                    do i=1,Hd_EBNR(n).order-1
                        if(usevar(i)==2) then!Carry
                            usevar(i)=0
                            usevar(i+1)=usevar(i+1)+1
                        end if
                    end do
                end do
            end if
        end do
        Hd_HdEC=HdECtemp
        do j=1,Hd_NState!Clean up
            do i=j,Hd_NState
                deallocate(HdECtemp(i,j).Array)
            end do
        end do
    end subroutine OriginShift

    !Load Hd expansion coefficient from file to HdEC
    !Optional: FileName: (default = 'Hd.CheckPoint') name of the input file
    subroutine ReadHdExpansionCoefficients(HdEC,FileName)
        type(d2PArray),dimension(Hd_NState,Hd_NState),intent(inout)::HdEC
        character*32,optional,intent(in)::FileName
        character*2::char2temp; character*28::char28temp
        integer::NState,NBasis,NOrder!The old Hd is not necessarily fitted under same condition
        integer::istate,jstate,i,j,order,location
        integer,allocatable,dimension(:)::indice
        real*8::dbletemp
        if(present(FileName)) then; open(unit=99,file=FileName,status='old')
        else; open(unit=99,file='Hd.CheckPoint',status='old'); end if
            read(99,'(A28,I2)')char28temp,NState!Get old Hd fitting condition
            read(99,*); read(99,*)dbletemp
            read(99,'(I5)')NOrder; allocate(indice(NOrder))
            NBasis=2!Get number of basis functions
            do
                read(99,'(A2)')char2temp; if(char2temp=='Hd') exit
                NBasis=NBasis+1
            end do
            NBasis=NBasis/2
            rewind 99!Number of basis functions has been gotten
            read(99,*)!Read old Hd
            do istate=1,NState
                do jstate=istate,NState
                    read(99,*)
                    do i=1,NBasis
                        read(99,*)dbletemp
                        read(99,'(I5)',advance='no')order
                        if(order>0) then
                            do j=1,order-1; read(99,'(I5)',advance='no')indice(j); end do; read(99,'(I5)')indice(order)
                        else
                            read(99,*)
                        end if
                        location=WhichExpansionBasis(order,indice(1:order))
                        if(location>0) then
                            HdEC(jstate,istate).Array(location)=dbletemp
                        else
                            write(*,*)'Warning: an old basis function does not belong to current basis space:'
                            do j=1,order-1; write(*,'(I5)',advance='no')indice(j); end do; write(*,'(I5)')indice(order)
                            write(*,*)'Its coefficient = ',dbletemp
                        end if
                    end do
                end do
            end do
        close(99)
    end subroutine ReadHdExpansionCoefficients

    !Write Hd expansion coefficient and expansion basis specification to file
    !Optional: FileName: (default = 'Hd.CheckPoint') name of the output file
    subroutine WriteHdExpansionCoefficients(HdEC,FileName)
        type(d2PArray),dimension(Hd_NState,Hd_NState),intent(in)::HdEC
        character*32,optional,intent(in)::FileName
        integer::istate,jstate,i,j
        if(present(FileName)) then; open(unit=99,file=FileName,status='replace')
        else; open(unit=99,file='Hd.CheckPoint',status='replace'); end if
    		write(99,'(A28,I2)')'Number of electronic states:',Hd_NState
            do istate=1,Hd_NState
                do jstate=istate,Hd_NState
                    write(99,'(A2,I2,I2)')'Hd',jstate,istate
                    do i=1,NHdExpansionBasis
                        write(99,*)HdEC(jstate,istate).Array(i)
                        if(Hd_EBNR(i).order>0) then
                            write(99,'(I5)',advance='no')Hd_EBNR(i).order
                            do j=1,Hd_EBNR(i).order-1; write(99,'(I5)',advance='no')Hd_EBNR(i).indice(j); end do; write(99,'(I5)')Hd_EBNR(i).indice(Hd_EBNR(i).order)
                        else
                            write(99,'(I5)')Hd_EBNR(i).order
                        end if
                    end do
                end do
            end do
        close(99)
    end subroutine WriteHdExpansionCoefficients

    !The value of n-th expansion basis function at some coordinate q
    real*8 function ExpansionBasis(q,n)
        real*8,dimension(Hd_intdim),intent(in)::q
        integer,intent(in)::n
        integer::i
        ExpansionBasis=1d0
        do i=1,Hd_EBNR(n).order
            ExpansionBasis=ExpansionBasis*q(Hd_EBNR(n).indice(i))
        end do
    end function ExpansionBasis

    !The value of ▽(n-th expansion basis function) at some coordinate q
    function ExpansionBasisGradient(q,n)
        real*8,dimension(Hd_intdim)::ExpansionBasisGradient
        real*8,dimension(Hd_intdim),intent(in)::q
        integer,intent(in)::n
        integer::m,i,OrderCount
        do m=1,Hd_intdim
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
        real*8,dimension(Hd_intdim,Hd_intdim)::ExpansionBasisHessian
        real*8,dimension(Hd_intdim),intent(in)::q
        integer,intent(in)::n
        integer::m1,m2,i,OrderCount1,OrderCount2
        do m1=1,Hd_intdim
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
            do m2=m1+1,Hd_intdim!Off-diagonal
                OrderCount1=0; OrderCount2=0
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

!------------- Fitting support --------------
	!To determine Hd, we need:
	!    Transformation of Hd expansion coefficient between HdEC form and vector c form
	!    Gradient of operators over c

    subroutine c2HdEC(c,HdEC,N)
        integer,intent(in)::N
        real*8,dimension(N),intent(in)::c
        type(d2PArray),dimension(Hd_NState,Hd_NState),intent(inout)::HdEC
        integer::istate,jstate,indice
        indice=1
        do istate=1,Hd_NState
            do jstate=istate,Hd_NState
                Hd_HdEC(jstate,istate).Array=c(indice:indice+NHdExpansionBasis-1)
                indice=indice+NHdExpansionBasis
            end do
        end do
    end subroutine c2HdEC
    subroutine HdEC2c(HdEC,c,N)
        integer,intent(in)::N
        type(d2PArray),dimension(Hd_NState,Hd_NState),intent(in)::HdEC
        real*8,dimension(N),intent(out)::c
        integer::istate,jstate,indice
        indice=1
        do istate=1,Hd_NState
            do jstate=istate,Hd_NState
                c(indice:indice+NHdExpansionBasis-1)=Hd_HdEC(jstate,istate).Array
                indice=indice+NHdExpansionBasis
            end do
        end do
	end subroutine HdEC2c
	
	!The value of ▽_cHd in diabatic representation at some coordinate q, where c is the expansion coefficient vector
    !f stores expansion basis function values at this q
    function dcHd_ByKnownf(f)
        real*8,dimension(NHdExpansionCoefficients,Hd_NState,Hd_NState)::dcHd_ByKnownf
        real*8,dimension(NHdExpansionBasis),intent(in)::f
        integer::i,j,indice
        indice=1
        do j=1,Hd_NState
            do i=j,Hd_NState
                dcHd_ByKnownf(1:indice-1,i,j)=0d0
                dcHd_ByKnownf(indice:indice+NHdExpansionBasis-1,i,j)=f
                indice=indice+NHdExpansionBasis
                dcHd_ByKnownf(indice:NHdExpansionCoefficients,i,j)=0d0
            end do
        end do
    end function dcHd_ByKnownf

    !The value of ▽_c▽Hd in diabatic representation at some coordinate q, where c is the expansion coefficient vector
    !fdT(i,:) stores the gradient of i-th expansion basis function
    function dcdHd_ByKnownfdT(fdT)
        real*8,dimension(NHdExpansionCoefficients,Hd_intdim,Hd_NState,Hd_NState)::dcdHd_ByKnownfdT
        real*8,dimension(NHdExpansionBasis,Hd_intdim),intent(in)::fdT
        integer::i,j,indice
        indice=1
        do j=1,Hd_NState
            do i=j,Hd_NState
                dcdHd_ByKnownfdT(1:indice-1,:,i,j)=0d0
                dcdHd_ByKnownfdT(indice:indice+NHdExpansionBasis-1,:,i,j)=fdT
                indice=indice+NHdExpansionBasis
                dcdHd_ByKnownfdT(indice:NHdExpansionCoefficients,:,i,j)=0d0
            end do
        end do
    end function dcdHd_ByKnownfdT

    !The value of ▽_cA in diabatic representation at some coordinate q, where c is the expansion coefficient vector
    !For A adopted here, i.e. (▽H)^2, we can use known ▽_c▽Hd & ▽Hd to calculate
    function dcAd_ByKnown(dcdHd,dHd)
        real*8,dimension(NHdExpansionCoefficients,Hd_NState,Hd_NState)::dcAd_ByKnown
        real*8,dimension(NHdExpansionCoefficients,Hd_intdim,Hd_NState,Hd_NState),intent(in)::dcdHd
        real*8,dimension(Hd_intdim,Hd_NState,Hd_NState),intent(in)::dHd
        dcAd_ByKnown=sy4matdotmulsy3(dcdHd,dHd,NHdExpansionCoefficients,Hd_intdim,Hd_NState)
        dcAd_ByKnown=dcAd_ByKnown+transpose3(dcAd_ByKnown,NHdExpansionCoefficients,Hd_NState,Hd_NState)
    end function dcAd_ByKnown
!------------------- End --------------------

!------------ Diabatic quantity -------------
    function Hd(q)!Return the value of Hd in diabatic representation at some coordinate q
        real*8,dimension(Hd_NState,Hd_NState)::Hd
        real*8,dimension(Hd_intdim),intent(in)::q
        integer::i,j
        real*8,dimension(NHdExpansionBasis)::f
        do i=1,NHdExpansionBasis
            f(i)=ExpansionBasis(q,i)
        end do
        forall(i=1:Hd_NState,j=1:Hd_NState,i>=j)
            Hd(i,j)=dot_product(f,Hd_HdEC(i,j).Array)
        end forall
    end function Hd

    function dHd(q)!Return the value of ▽Hd in diabatic representation at some coordinate q
        real*8,dimension(Hd_intdim,Hd_NState,Hd_NState)::dHd
        real*8,dimension(Hd_intdim),intent(in)::q
        integer::i,j
        real*8,dimension(Hd_intdim,NHdExpansionBasis)::fd
        do i=1,NHdExpansionBasis
            fd(:,i)=ExpansionBasisGradient(q,i)
        end do
        forall(i=1:Hd_NState,j=1:Hd_NState,i>=j)
            dHd(:,i,j)=matmul(fd,Hd_HdEC(i,j).Array)
        end forall
    end function dHd

    function ddHd(q)!Return the value of ▽▽Hd in diabatic representation at some coordinate q
        real*8,dimension(Hd_intdim,Hd_intdim,Hd_NState,Hd_NState)::ddHd
        real*8,dimension(Hd_intdim),intent(in)::q
        integer::i,j,k
        real*8,dimension(Hd_intdim,Hd_intdim,NHdExpansionBasis)::fdd
        do i=1,NHdExpansionBasis
            fdd(:,:,i)=ExpansionBasisHessian(q,i)
        end do
        forall(i=1:Hd_NState,j=1:Hd_NState,k=1:Hd_intdim,i>=j)
            ddHd(:,k,i,j)=matmul(fdd(:,k,:),Hd_HdEC(i,j).Array)
        end forall
    end function ddHd

    !The value of Hd in diabatic representation and expansion basis functions at some coordinate q
    subroutine Hd_f(Hd,f,q)
        real*8,dimension(Hd_NState,Hd_NState),intent(out)::Hd
        real*8,dimension(NHdExpansionBasis),intent(out)::f
        real*8,dimension(Hd_intdim),intent(in)::q
        integer::i,j
        do i=1,NHdExpansionBasis
            f(i)=ExpansionBasis(q,i)
        end do
        forall(i=1:Hd_NState,j=1:Hd_NState,i>=j)
            Hd(i,j)=dot_product(f,Hd_HdEC(i,j).Array)
        end forall
    end subroutine Hd_f

    !The value of ▽Hd in diabatic representation and expansion basis function gradient at some coordinate q
    !fd(:,i) stores the gradient of i-th expansion basis function
    subroutine dHd_fd(dHd,fd,q)
        real*8,dimension(Hd_intdim,Hd_NState,Hd_NState),intent(out)::dHd
        real*8,dimension(Hd_intdim,NHdExpansionBasis),intent(out)::fd
        real*8,dimension(Hd_intdim),intent(in)::q
        integer::i,j
        do i=1,NHdExpansionBasis
            fd(:,i)=ExpansionBasisGradient(q,i)
        end do
        forall(i=1:Hd_NState,j=1:Hd_NState,i>=j)
            dHd(:,i,j)=matmul(fd,Hd_HdEC(i,j).Array)
        end forall
    end subroutine dHd_fd
!------------------- End --------------------

!------------ Adiabatic quantity ------------
    !Compute adiabatic quantity from Hd at some coordinate q

    function AdiabaticEnergy(q)!Return adiabatic energy
        real*8,dimension(Hd_NState)::AdiabaticEnergy
        real*8,dimension(Hd_intdim),intent(in)::q
        real*8,dimension(Hd_NState,Hd_NState)::phi
        phi=Hd(q)
        call My_dsyev('N',phi,AdiabaticEnergy,Hd_NState)
    end function AdiabaticEnergy

    function AdiabaticdH(q)!Return adiabatic gradient
        real*8,dimension(Hd_intdim,Hd_NState,Hd_NState)::AdiabaticdH
        real*8,dimension(Hd_intdim),intent(in)::q
        real*8,dimension(Hd_NState)::energy
        real*8,dimension(Hd_NState,Hd_NState)::phi
        phi=Hd(q)
        call My_dsyev('V',phi,energy,Hd_NState)
        AdiabaticdH=dHd(q)
        AdiabaticdH=sy3UnitaryTransformation(AdiabaticdH,phi,Hd_intdim,Hd_NState)
    end function AdiabaticdH

    function AdiabaticddH(q)!Return adiabatic Hessian
        real*8,dimension(Hd_intdim,Hd_intdim,Hd_NState,Hd_NState)::AdiabaticddH
        real*8,dimension(Hd_intdim),intent(in)::q
        real*8,dimension(Hd_intdim,Hd_NState,Hd_NState)::dH,M
        real*8,dimension(Hd_NState)::energy
        real*8,dimension(Hd_NState,Hd_NState)::phi
        phi=Hd(q)
        call My_dsyev('V',phi,energy,Hd_NState)
        dH=dHd(q)
        dH=sy3UnitaryTransformation(dH,phi,Hd_intdim,Hd_NState)
        M=deigvec_ByKnowneigval_dA(energy,dH,Hd_intdim,Hd_NState)
        AdiabaticddH=asy3matdirectmulsy3(M,dH,Hd_intdim,Hd_intdim,Hd_NState)
        AdiabaticddH=sy4UnitaryTransformation(ddHd(q),phi,Hd_intdim,Hd_intdim,Hd_NState)&
            -AdiabaticddH-transpose4(AdiabaticddH,Hd_intdim,Hd_intdim,Hd_NState,Hd_NState)
    end function AdiabaticddH

    !energy harvests adiabatic energy, dH harvests ▽H_a
    subroutine AdiabaticEnergy_dH(q,energy,dH)
        real*8,dimension(Hd_intdim),intent(in)::q
        real*8,dimension(Hd_NState),intent(out)::energy
        real*8,dimension(Hd_intdim,Hd_NState,Hd_NState),intent(out)::dH
        real*8,dimension(Hd_NState,Hd_NState)::phi
        phi=Hd(q)
        call My_dsyev('V',phi,energy,Hd_NState)
        dH=dHd(q)
        dH=sy3UnitaryTransformation(dH,phi,Hd_intdim,Hd_NState)
    end subroutine AdiabaticEnergy_dH

    !phi harvests adiabatic states in diabatic representation
    subroutine AdiabaticEnergy_dH_State(q,energy,dH,phi)
        real*8,dimension(Hd_intdim),intent(in)::q
        real*8,dimension(Hd_NState),intent(out)::energy
        real*8,dimension(Hd_intdim,Hd_NState,Hd_NState),intent(out)::dH
        real*8,dimension(Hd_NState,Hd_NState),intent(out)::phi
        phi=Hd(q)
        call My_dsyev('V',phi,energy,Hd_NState)
        dH=dHd(q)
        dH=sy3UnitaryTransformation(dH,phi,Hd_intdim,Hd_NState)
    end subroutine AdiabaticEnergy_dH_State

    !f harvests expansion basis function values 
    subroutine AdiabaticEnergy_State_f(q,energy,phi,f)
        real*8,dimension(Hd_intdim),intent(in)::q
        real*8,dimension(Hd_NState),intent(out)::energy
        real*8,dimension(Hd_NState,Hd_NState),intent(out)::phi
        real*8,dimension(NHdExpansionBasis),intent(out)::f
        call Hd_f(phi,f,q)
        call My_dsyev('V',phi,energy,Hd_NState)
    end subroutine AdiabaticEnergy_State_f

    !fd harvests expansion basis function gradient values, fd(:,i) = the gradient of i-th expansion basis function
    subroutine AdiabaticEnergy_dH_State_f_fd(q,energy,dH,phi,f,fd)
        real*8,dimension(Hd_intdim),intent(in)::q
        real*8,dimension(Hd_NState),intent(out)::energy
        real*8,dimension(Hd_intdim,Hd_NState,Hd_NState),intent(out)::dH
        real*8,dimension(Hd_NState,Hd_NState),intent(out)::phi
        real*8,dimension(NHdExpansionBasis),intent(out)::f
        real*8,dimension(Hd_intdim,NHdExpansionBasis),intent(out)::fd
        call Hd_f(phi,f,q)
        call My_dsyev('V',phi,energy,Hd_NState)
        call dHd_fd(dH,fd,q)
        dH=sy3UnitaryTransformation(dH,phi,Hd_intdim,Hd_NState)
    end subroutine AdiabaticEnergy_dH_State_f_fd
!------------------- End --------------------

!---------- Nondegenerate quantity ----------
    !When Hamiltonian is almost degenerate, adiabatic basis suffers from singularity
    !This causes no trouble to H, but operators other than H suffer a lot (e.g. ▽H)
    !Diabatz is always good, but it cannot be determined without knowledge of neighbourhood,
    !sometimes we ask for a nondegenerate representation based only on single point information
    !Choose (▽H)^2 representation: it is good enough for fitting a C1 molecule when not exploring too much seam space
    !    (▽H)^2 is Hermitian so there exists exactly one (▽H)^2 representation
    !    It reduces to s orthogonal to h for 2-fold degeneracy case, where s is the average force of 2 PESs
    !Note it is only appropriate around conical intersection, because || ▽H || -> 0 at asymptote

    !Compute nondegenerate representation quantity from Hd at some coordinate q
    !Optional: DegenerateThreshold: if present, check whether some eigenvalues closer than DegenerateThreshold

    !dim x NState x NState 3-order tensor ▽H
    !Transform ▽H to (▽H)^2 representation
    !eigval harvests eigen values of (▽H)^2
	!eigvec harvests eigen vectors of (▽H)^2 in representation same to input ▽H
    subroutine NondegenerateRepresentation(dH,eigval,eigvec,dim,NState,DegenerateThreshold)
        integer,intent(in)::dim,NState
        real*8,dimension(dim,NState,NState),intent(inout)::dH
        real*8,dimension(NState),intent(out)::eigval
		real*8,dimension(NState,NState),intent(out)::eigvec
		real*8,intent(in),optional::DegenerateThreshold
        logical::degenerate
        eigvec=sy3matdotmul(dH,dH,dim,NState); call My_dsyev('V',eigvec,eigval,NState)
		dH=sy3UnitaryTransformation(dH,eigvec,dim,NState)
		if(present(DegenerateThreshold)) then
            call CheckDegeneracy(degenerate,DegenerateThreshold,eigval,NState)
			if(degenerate) write(*,*)'Warning: nondegenerate representation has eigenvalues closer than',DegenerateThreshold
		end if
    end subroutine NondegenerateRepresentation

    !H harvests H_nd, dH harvests ▽H_nd
    subroutine NondegenerateH_dH(q,H,dH,DegenerateThreshold)
        real*8,dimension(Hd_intdim),intent(in)::q
        real*8,dimension(Hd_NState,Hd_NState),intent(out)::H
        real*8,dimension(Hd_intdim,Hd_NState,Hd_NState),intent(out)::dH
        real*8,intent(in),optional::DegenerateThreshold
        real*8,dimension(Hd_NState)::eigval
        real*8,dimension(Hd_NState,Hd_NState)::phi
        H=Hd(q); dH=dHd(q)
        if(present(DegenerateThreshold)) then
            call NondegenerateRepresentation(dH,eigval,phi,Hd_intdim,Hd_NState,DegenerateThreshold)
        else
            call NondegenerateRepresentation(dH,eigval,phi,Hd_intdim,Hd_NState)
        end if
        H=matmul(transpose(phi),matmul(H,phi))
    end subroutine NondegenerateH_dH

    !phi harvests the nondegenerate states in diabatic representation
    !f harvests expansion basis function values
    !fd harvests expansion basis function gradient values, fd(:,i) = the gradient of i-th expansion basis function
    subroutine NondegenerateH_dH_State_f_fd(q,H,dH,phi,f,fd,DegenerateThreshold)
        real*8,dimension(Hd_intdim),intent(in)::q
        real*8,dimension(Hd_NState,Hd_NState),intent(out)::H
        real*8,dimension(Hd_intdim,Hd_NState,Hd_NState),intent(out)::dH
        real*8,dimension(Hd_NState,Hd_NState),intent(out)::phi
        real*8,dimension(NHdExpansionBasis),intent(out)::f
        real*8,dimension(Hd_intdim,NHdExpansionBasis),intent(out)::fd
        real*8,intent(in),optional::DegenerateThreshold
        real*8,dimension(Hd_NState)::eigval
        call Hd_f(H,f,q); call dHd_fd(dH,fd,q)
        if(present(DegenerateThreshold)) then
            call NondegenerateRepresentation(dH,eigval,phi,Hd_intdim,Hd_NState,DegenerateThreshold)
        else
            call NondegenerateRepresentation(dH,eigval,phi,Hd_intdim,Hd_NState)
        end if
        H=matmul(transpose(phi),matmul(H,phi))
    end subroutine NondegenerateH_dH_State_f_fd

    !eigval harvests the eigenvalues corresponding to phi
    !dHd harvests ▽H_d
    subroutine NondegenerateH_dH_eigval_State_dHd_f_fd(q,H,dH,eigval,phi,dHd,f,fd,DegenerateThreshold)
        real*8,dimension(Hd_intdim),intent(in)::q
        real*8,dimension(Hd_NState,Hd_NState),intent(out)::H
        real*8,dimension(Hd_intdim,Hd_NState,Hd_NState),intent(out)::dH
        real*8,dimension(Hd_NState),intent(out)::eigval
        real*8,dimension(Hd_NState,Hd_NState),intent(out)::phi
        real*8,dimension(Hd_intdim,Hd_NState,Hd_NState),intent(out)::dHd
        real*8,dimension(NHdExpansionBasis),intent(out)::f
        real*8,dimension(Hd_intdim,NHdExpansionBasis),intent(out)::fd
        real*8,intent(in),optional::DegenerateThreshold
        call Hd_f(H,f,q); call dHd_fd(dHd,fd,q); dH=dHd
        if(present(DegenerateThreshold)) then
            call NondegenerateRepresentation(dH,eigval,phi,Hd_intdim,Hd_NState,DegenerateThreshold)
        else
            call NondegenerateRepresentation(dH,eigval,phi,Hd_intdim,Hd_NState)
        end if
        H=matmul(transpose(phi),matmul(H,phi))
    end subroutine NondegenerateH_dH_eigval_State_dHd_f_fd
!------------------- End --------------------

end module DiabaticHamiltonian