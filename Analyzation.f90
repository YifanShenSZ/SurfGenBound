!Provide analyzation of the fitted potential energy surface:
!    Compute H & ▽H in diabatic & adiabatic & nondegenerate representation on specified geometries
!    Minimum search and vibration analysis of specified adiabatic state
!    Mex search and gh orthogonalization along conical intersection seam between specified adiabatic states
module Analyzation
    use Basic
    use DiabaticHamiltonian
    implicit none

!Analyzation module only variable
    character*32::AnalyzeJob
    integer::InterestingState
    
contains
!Read the input file for Analyzation: AnalyzeInput 
subroutine ReadAnalyzeInput()
    open(unit=99)
end subroutine ReadAnalyzeInput

subroutine Analyze()
end subroutine Analyze

subroutine MinimumSearch()
end subroutine MinimumSearch

subroutine MexSearch()
end subroutine MexSearch

subroutine Evaluate()
end subroutine Evaluate

!---------- Auxiliary routine ----------
!----------------- End -----------------

end module Analyzation