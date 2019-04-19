!Generate NadVibS input files
module NadVibS
    use Basic
    use ESSInput
    use DiabaticHamiltonian
    use Analyze
    implicit none

contains
subroutine GenerateNadVibsInput()
    integer::ip,istate,jstate,iorder
    real*8,dimension(InternalDimension)::qPrecursor,qSuccessor,freqrPrecursor,freqrSuccessor,freqiPrecursor,freqiSuccessor
    real*8,dimension(InternalDimension,InternalDimension)::HPrecursor,HSuccessor
    real*8,dimension(InternalDimension,CartesianDimension)::BPrecursor,BSuccessor
    real*8,dimension(InternalDimension,InternalDimension,NStates,NStates)::Htemp
    type(Data),allocatable,dimension(:)::pointtemp
    !Definition of dshift and Tshift see Schuurman & Yarkony 2008 JCP 128 eq. (12)
    real*8,dimension(InternalDimension)::dshift
    real*8,dimension(InternalDimension,InternalDimension)::Tshift
    !Precursor
    call WilsonBMatrixAndInternalCoordinateq(BPrecursor,qPrecursor,reshape(MoleculeDetail.RefConfig,[CartesianDimension]),InternalDimension,CartesianDimension)
    call ReadESSHessian(HPrecursor,InternalDimension)
    call VibrationAnalysis(freqrPrecursor,freqiPrecursor,HPrecursor,InternalDimension,BPrecursor,CartesianDimension,MoleculeDetail.mass,NAtoms)
    if(maxval(abs(freqiPrecursor))>1d-14) write(*,*)'Warning: imaginary frequency found for precursor'
    !Successor
        !Allocate storage space
            allocate(pointtemp(NPoints))
            do ip=1,NPoints
                allocate(pointtemp(ip).geom(CartesianDimension))
                allocate(pointtemp(ip).energy(NStates))
                allocate(pointtemp(ip).dH(CartesianDimension,NStates,NStates))
            end do
    call ReadESSData(pointtemp,NPoints)
    call WilsonBMatrixAndInternalCoordinateq(BSuccessor,qSuccessor,pointtemp(IndexReference).geom,InternalDimension,CartesianDimension)
!This version directly use Hessian at reference geometry
    Htemp=AdiabaticddH(qSuccessor)
    HSuccessor=Htemp(:,:,1,1)
!I will write a version shift the reference to ground state minimum of Hd someday
    call VibrationAnalysis(freqrSuccessor,freqiSuccessor,HSuccessor,InternalDimension,BSuccessor,CartesianDimension,MoleculeDetail.mass,NAtoms)
    if(maxval(abs(freqiSuccessor))>1d-14) write(*,*)'Warning: imaginary frequency found for successor'
    dshift=matmul(transpose(HPrecursor),qSuccessor-qPrecursor)
    Tshift=matmul(transpose(HPrecursor),HSuccessor)
    open(unit=99,file='nadvibs.in',status='replace')
        write(99,'(A59)')'Angular frequency in atomic unit of each vibrational basis:'
        write(99,*)freqrSuccessor
        do istate=1,NStates
            do jstate=istate,NStates
                do iorder=0,NOrder
                    write(99,'(A2,I2,I2,A15,I2)')'Hd',jstate,istate,'Expansion order',iorder
                    write(99,*)HdEC(jstate,istate).Order(iorder).Array
                end do
            end do
        end do
        write(99,'(A63)')'Angular frequency in atomic unit of each precursor normal mode:'
        write(99,*)freqrPrecursor
        write(99,'(A13)')'Shift Vector:'
        write(99,*)dshift
        write(99,'(A22)')'Transformation Matrix:'
        write(99,*)Tshift
    close(99)
end subroutine GenerateNadVibsInput

end module NadVibs