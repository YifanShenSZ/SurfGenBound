!Read concatenated electronic structure software (ESS) output
!In this specific code, ESS = Columbus7
module ESSInput
    use Basic
    implicit none

!Parameter
    character*32,parameter::ElectronicStructureSoftware='Columbus7'

contains
subroutine ReadESSOutput(point,NPoints)
    integer,intent(in)::NPoints
    type(Data),dimension(NPoints),intent(inout)::point
    character*1::CharTemp1,CharTemp2
    character*128::source
    integer::ip,istate,jstate
    write(*,*)'Electronic structure software = '//ElectronicStructureSoftware
    source='geom.all'
    call ColumbusGeometry(source,point,NPoints)
    source='energy.all'
    call ColumbusEnergy(source,point,NPoints)
    do istate=1,NStates
        write(CharTemp1,'(I1)')istate
        source='cartgrd.drt1.state'//CharTemp1//'.all'
        call ColumbusGradient(source,istate,istate,point,NPoints)
        do jstate=istate+1,NStates
            write(CharTemp2,'(I1)')jstate
            source='cartgrd.nad.drt1.state'//CharTemp1//'.drt1.state'//CharTemp2//'.all'
            call ColumbusGradient(source,jstate,istate,point,NPoints)
        end do
    end do
end subroutine ReadESSOutput

subroutine DefineInternalCoordinate()
    write(*,*)'Define internal coordinate according to '//ElectronicStructureSoftware
    call Columbusintcfl(InternalDimension)
end subroutine DefineInternalCoordinate

subroutine ColumbusGeometry(source,point,NPoints)
    character*128,intent(in)::source
    integer,intent(in)::NPoints
    type(Data),dimension(NPoints),intent(inout)::point
    integer::ip,iatm
    character*2::CharTemp
    real*8::DbTemp1,DbTemp2
    open(unit=99,file=source,status='old')
        do ip=1,NPoints
            do iatm=1,NAtoms
                read(99,'(1x,a2,2x,f5.1,4f14.8)')CharTemp,DbTemp1,point(ip).geom(3*iatm-2:3*iatm),DbTemp2
            end do
        end do
    close(99)
end subroutine ColumbusGeometry

subroutine ColumbusEnergy(source,point,NPoints)
    character*128,intent(in)::source
    integer,intent(in)::NPoints
    type(Data),dimension(NPoints),intent(inout)::point
    integer::ip
    open(unit=99,file=source,status='old')
        do ip=1,NPoints
            read(99,'(f14.8,1x,f14.8)')point(ip).energy
        end do
    close(99)
end subroutine ColumbusEnergy

subroutine ColumbusGradient(source,state1,state2,point,NPoints)
    character*128,intent(in)::source
    integer,intent(in)::state1,state2,NPoints
    type(Data),dimension(NPoints),intent(inout)::point
    integer::ip,iatm
    open(unit=99,file=source,status='old')
        do ip=1,NPoints
            do iatm=1,NAtoms
                read(99,'(3f15.8)')point(ip).dH(3*iatm-2:3*iatm,state1,state2)
            end do
        end do
    close(99)
end subroutine ColumbusGradient

subroutine Columbusintcfl(intdim)
    integer,intent(out)::intdim!Return the dimension of internal space
    character*10,allocatable,dimension(:)::MotionType
    character*24::CharTemp24
    integer::i,j,k,NLines
    integer,allocatable,dimension(:)::KLine
    real*8::DbTemp
    open(unit=99,file='intcfl',status='old')
        !Get how many lines are the definition for internal coordinates & how many internal coordinates there are
            NLines=0
            intdim=0
            read(99,*)!First line is always 'TEXAS'
            do while(.true.)
                read(99,'(A24)')CharTemp24
                if (index(CharTemp24,'STRE')==0.and.index(CharTemp24,'BEND')==0.and.index(CharTemp24,'TORS')==0) exit
                NLines=NLines+1
                if(scan(CharTemp24,'K')==1) intdim=intdim+1
            end do
            rewind 99
        !Get whether a line is the start of a new internal coordinate, and what type of motion a line stands for
            allocate(KLine(intdim+1))
            KLine(intdim+1)=NLines+1
            allocate(MotionType(NLines))
            j=1
            read(99,*)!First line is always 'TEXAS'
            do i=1,NLines
                read(99,'(A24)')CharTemp24
                if(scan(CharTemp24,'K')==1) then
                    KLine(j)=i
                    j=j+1
                end if
                if(index(CharTemp24,'STRE')>0) then
                    MotionType(i)='stretching'
                else if(index(CharTemp24,'BEND')>0) then
                    MotionType(i)='bending'
                else if(index(CharTemp24,'TORS')>0) then
                    MotionType(i)='torsion'
                end if
            end do
            rewind 99
        !Finally read Columbus internal coordinate definition
            allocate(GeometryTransformation_IntCDef(intdim))
            k=1!Counter for line
            read(99,*)!First line is always 'TEXAS'
            do i=1,intdim
                GeometryTransformation_IntCDef(i).NMotions=KLine(i+1)-KLine(i)
                allocate(GeometryTransformation_IntCDef(i).motion(GeometryTransformation_IntCDef(i).NMotions))
                if(GeometryTransformation_IntCDef(i).NMotions==1) then
                    GeometryTransformation_IntCDef(i).motion(1).type=MotionType(k)
                    GeometryTransformation_IntCDef(i).motion(1).coeff=1d0!Only 1 motion
                    select case(MotionType(k))
                        case('stretching')
                            allocate(GeometryTransformation_IntCDef(i).motion(1).atom(2))
                            read(99,'(A28,I5,1x,I9)')CharTemp24,GeometryTransformation_IntCDef(i).motion(1).atom
                        case('bending')
                                allocate(GeometryTransformation_IntCDef(i).motion(1).atom(3))
                                read(99,'(A28,I5,1x,I9,1x,I9)')CharTemp24,GeometryTransformation_IntCDef(i).motion(1).atom(1),&
                                    GeometryTransformation_IntCDef(i).motion(1).atom(3),GeometryTransformation_IntCDef(i).motion(1).atom(2)
                        case('torsion')
                            allocate(GeometryTransformation_IntCDef(i).motion(1).atom(4))
                            read(99,'(A28,I5,1x,I9,1x,I9,1x,I9)')CharTemp24,GeometryTransformation_IntCDef(i).motion(1).atom
                        case default!Throw a warning
                            write(*,'(1x,A51,1x,A10)')'Program abort: unsupported internal coordinate type',MotionType(k)
                            stop
                    end select
                    k=k+1
                else
                    DbTemp=0d0
                    do j=1,GeometryTransformation_IntCDef(i).NMotions
                        GeometryTransformation_IntCDef(i).motion(j).type=MotionType(k)
                        select case(MotionType(k))
                            case('stretching')
                                allocate(GeometryTransformation_IntCDef(i).motion(j).atom(2))
                                read(99,'(A10,F10.7,8x,I6,1x,I9)')CharTemp24,&
                                    GeometryTransformation_IntCDef(i).motion(j).coeff,GeometryTransformation_IntCDef(i).motion(j).atom
                            case('bending')
                                allocate(GeometryTransformation_IntCDef(i).motion(j).atom(3))
                                read(99,'(A10,F10.7,8x,I6,1x,I9,1x,I9)')CharTemp24,GeometryTransformation_IntCDef(i).motion(j).coeff,&
                                    GeometryTransformation_IntCDef(i).motion(j).atom(1),GeometryTransformation_IntCDef(i).motion(j).atom(3),GeometryTransformation_IntCDef(i).motion(j).atom(2)
                            case('torsion')
                                allocate(GeometryTransformation_IntCDef(i).motion(j).atom(4))
                                read(99,'(A10,F10.7,8x,I6,1x,I9,1x,I9,1x,I9)')CharTemp24,&
                                    GeometryTransformation_IntCDef(i).motion(j).coeff,GeometryTransformation_IntCDef(i).motion(j).atom
                            case default!Throw a warning
                                write(*,'(1x,A51,1x,A10)')'Program abort: unsupported internal coordinate type',MotionType(k)
                                stop
                        end select
                        k=k+1
                        DbTemp=DbTemp+GeometryTransformation_IntCDef(i).motion(j).coeff*GeometryTransformation_IntCDef(i).motion(j).coeff
                    end do
                    DbTemp=Sqrt(DbTemp)
                    forall(j=1:GeometryTransformation_IntCDef(i).NMotions)
                        GeometryTransformation_IntCDef(i).motion(j).coeff=GeometryTransformation_IntCDef(i).motion(j).coeff/DbTemp
                    end forall
                end if
            end do
    close(99)
    !Clean up
        deallocate(MotionType)
        deallocate(KLine)
end subroutine Columbusintcfl

end module ESSInput