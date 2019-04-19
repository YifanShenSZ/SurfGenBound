!Generate NadVibS input files
module NadVibS
    use Basic
    use ESSInput
    use DiabaticHamiltonian
    implicit none

contains
subroutine GenerateNadVibsInput()
    integer::ip,istate,jstate,iorder
    real*8,dimension(InternalDimension)::qPrecursor,qSuccesor
    real*8,dimension(InternalDimension,CartesianDimension)::BPrecursor,BSuccesor
    type(Data),allocatable,dimension(:)::pointtemp
    !Definition of dshift and Tshift see Schuurman & Yarkony 2008 JCP 128 eq. (12)
    real*8,dimension(InternalDimension)::dshift
    real*8,dimension(InternalDimension,InternalDimension)::Tshift
    !Read successor reference Cartesian geometry. We already have precursor Cartesian geometry in MoleculeDetail.RefConfig
        !Allocate storage space
            allocate(pointtemp(NPoints))
            do ip=1,NPoints
                allocate(pointtemp(ip).geom(CartesianDimension))
                allocate(pointtemp(ip).energy(NStates))
                allocate(pointtemp(ip).dH(CartesianDimension,NStates,NStates))
            end do
            allocate(ReferencePoint.geom(CartesianDimension))
        call ReadESSOutput(pointtemp,NPoints)
        ReferencePoint.geom=pointtemp(IndexReference).geom
        deallocate(pointtemp)!Clean up
    !Allocate storage space
    call WilsonBMatrixAndInternalCoordinateq(BPrecursor,qPrecursor,reshape(MoleculeDetail.RefConfig,[CartesianDimension]),InternalDimension,CartesianDimension)
    call WilsonBMatrixAndInternalCoordinateq(BSuccesor,qSuccesor,ReferencePoint.geom,InternalDimension,CartesianDimension)
    open(unit=99,file='nadvibs.in',status='replace')
        write(99,'(A48)')'hbar * omega in Hatree of each vibrational basis'
        
        do istate=1,NStates
            do jstate=istate,NStates
                do iorder=0,NOrder
                    write(99,'(A2,I2,I2,A15,I2)')'Hd',jstate,istate,'Expansion order',iorder
                    write(99,*)HdEC(jstate,istate).Order(iorder).Array
                end do
            end do
        end do

    close(99)
end subroutine GenerateNadVibsInput

end module NadVibs