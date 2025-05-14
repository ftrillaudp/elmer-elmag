FUNCTION getElectricalResistivity(Model, n, arg) RESULT(rho)
  !!! Frederic Trillaud <ftrillaudp@gmail.com> - May 2025
  !!! elmerf90 -o electricalResistivity.so electricalResistivity.F90
  !!! Compute the local electric resistivity of the superconductor rho 

  ! Elmer module
  USE DefUtils

  IMPLICIT NONE
  TYPE(Model_t) :: model
  INTEGER :: n
  REAL(KIND=dp) :: arg(*)
  REAL(KIND=dp) :: Je, Jex, Jey, Jez, nv, Ec, Jc, rho
  LOGICAL :: gotIt, visu

  !!! Get parameters from sif file:
  ! variables needed inside function
  TYPE(ValueList_t), POINTER :: material
  ! get pointer on list for material
  material => GetMaterial()
  IF (.NOT. ASSOCIATED(material)) THEN
    CALL Fatal('getElectricalConductivity', 'No material found')
  END IF
  ! read in reference n-value
  nv = GetConstReal(material, 'N-Value', gotIt)
  IF (.NOT. gotIt) THEN
    CALL Fatal('getElectricalConductivity', 'N-Value')
  END IF
  Ec = GetConstReal(material, 'Critical Electrical Field', gotIt)
  IF (.NOT. gotIt) THEN
    CALL Fatal('getElectricalConductivity', 'Critical Electrical Field')
  END IF
  Jc = GetConstReal(material, 'Critical Current Density', gotIt)
  IF (.NOT. gotIt) THEN
    CALL Fatal('getElectricalConductivity', 'Critical Current Density')
  END IF
  
  visu = .FALSE.

  !!! Get the variables from the input:
  Jex = arg(1)
  Jey = arg(2)
  Jez = arg(3)
  Je = SQRT(Jex**2+Jey**2+Jez**2)
  
  !!! Electric resistivity:
  ! Superconductor
  rho = (Ec / Jc**nv)*(Je**(nv-1))
  
  IF (visu) THEN
    PRINT 1, rho
    1  FORMAT(' rho_sc: ', EN12.3)
  END IF

END FUNCTION getElectricalResistivity
