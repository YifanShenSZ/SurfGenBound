!Read concatenated electronic structure software output
!In this specific code, electronic structure software = Columbus7
module ElectronicStructure
    use Basic
    implicit none

!Parameter
    character*32,parameter::ElectronicStructureSoftware='Columbus7'

contains
!Top level routines, providing standard interface for other modules to call

subroutine ReadElectronicStructureData(point,NPoints)!Fill in point.geom & energy & dH
    integer,intent(in)::NPoints
    type(Data),dimension(NPoints),intent(inout)::point
    character*1::chartemp1,chartemp2
    character*128::source
    integer::ip,istate,jstate
    source='geom.all'
    call ColumbusGeometry(source,point,NPoints)
    source='energy.all'
    call ColumbusEnergy(source,point,NPoints)
    do istate=1,NState
        write(chartemp1,'(I1)')istate
        source='cartgrd.drt1.state'//chartemp1//'.all'
        call ColumbusGradient(source,istate,istate,point,NPoints)
        do jstate=istate+1,NState
            write(chartemp2,'(I1)')jstate
            source='cartgrd.nad.drt1.state'//chartemp1//'.drt1.state'//chartemp2//'.all'
            call ColumbusGradient(source,jstate,istate,point,NPoints)
        end do
    end do
end subroutine ReadElectronicStructureData

subroutine ReadElectronicStructureHessian(H,intdim)
    integer,intent(in)::intdim
    real*8,dimension(intdim,intdim),intent(out)::H
    call ColumbusHessian(H,intdim)
end subroutine ReadElectronicStructureHessian

!---------- Bottom level ----------
    !Actually dealing with certain electronic structure software

    subroutine ColumbusGeometry(source,point,NPoints)
        character*128,intent(in)::source
        integer,intent(in)::NPoints
        type(Data),dimension(NPoints),intent(inout)::point
        integer::ip,iatm
        character*2::chartemp
        real*8::dbletemp1,dbletemp2
        open(unit=99,file=source,status='old')
            do ip=1,NPoints
                do iatm=1,MoleculeDetail.NAtoms
                    read(99,*)chartemp,dbletemp1,point(ip).geom(3*iatm-2:3*iatm),dbletemp2
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
                read(99,*)point(ip).energy
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
                do iatm=1,MoleculeDetail.NAtoms
                    read(99,*)point(ip).dH(3*iatm-2:3*iatm,state1,state2)
                end do
            end do
        close(99)
    end subroutine ColumbusGradient
    
    !Cris provides me with Fmat rather than hessian,
    !so this routine reads Columbus Fmat file rather than hessian file
    subroutine ColumbusHessian(H,intdim)
        integer,intent(in)::intdim
        real*8,dimension(intdim,intdim),intent(out)::H
        integer::i
        open(unit=99,file='Fmat',status='old')
            do i=1,intdim
                read(99,*)H(i,1:i)
            end do
        close(99)
        !The internal coordinate and vibration routines of Columbus use weird unit:
        !    energy in 10^-18 J, length in A (to be continued)
        H=H/4.35974417d0! 1 Hatree = 4.35974417 * 10^-18 J
        do i=1,intdim
            if(GeometryTransformation_IntCoordDef(i).motion(1).type=='stretching') then
                H(:,i)=H(:,i)/AInAU
                H(i,:)=H(i,:)/AInAU
            end if
        end do
    end subroutine ColumbusHessian
!-------------- End ---------------

end module ElectronicStructure