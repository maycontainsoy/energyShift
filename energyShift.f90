PROGRAM rmsFortran

  IMPLICIT NONE

  !......................................................................!
  INTEGER :: i, j     ! Counters
  INTEGER :: nlines   ! Number of lines in file
  INTEGER :: discard  ! Number of rows to discard
  INTEGER :: n

  REAL :: sum1      ! <(Ephf - Esm)^2>
  REAL :: sum2      ! <(Ephf - Esm)>^2
  REAL :: total     ! sum1 - sum2 = sigma^2
  REAL :: totalOld  ! Old variance
  REAL :: totalNew  ! New variance
  REAL :: criteria  ! Criteria for discarding data
  REAL :: deltaE    ! Shift in absolute PHF energies
  REAL :: compare

  CHARACTER (LEN=30) :: fileName
  CHARACTER (LEN=30) :: outFile

  REAL, ALLOCATABLE :: energies(:,:)
  !......................................................................!

  PRINT *, '!......................................................................!'
  PRINT *, '!                   Program: rmsFortran_v0.f90                         !'
  PRINT *, '!                    Created: September 2020                           !'
  PRINT *, '!                  Created by: Stephanie Lauber                        !'
  PRINT *, '!                                                                      !'
  PRINT *, '! Add description here                                                 !'
  PRINT *, '!......................................................................!'

  ! Initialize variables
  nlines  = 0
  sum1    = 0.
  sum2    = 0.
  total   = 0.
  discard = 0

  PRINT *, 'Enter name of pre-processed file: '
  READ (*,*) fileName
  fileName = TRIM(fileName)

  PRINT *, 'Enter name for output file: '
  READ (*,*) outFile
  outFile = TRIM(outFile)

  PRINT *, 'Enter a value for energy shift (delta E): '
  READ (*,*) deltaE

  PRINT *, 'Enter n for n*sigma discard: '
  READ (*,*) n

  OPEN (1,file = fileName)
  ! Read in number of rows in .txt file
  DO
    READ (1,*,END=10)
    nlines = nlines + 1
  END DO
  10 CLOSE (1)

  PRINT *, 'Number of data pairs: ', nlines

  ! Allocate dimensions to energy array
  ! Col Data
  ! 1   PHF absolute energy
  ! 2   PHF excitation energy
  ! 3   FCI absolute energy
  ! 4   FCI excitation energy
  ! 5   Angular Momentum
  ALLOCATE(energies(nlines,5))

  OPEN (1,file = fileName)
  ! Read energies from file into array
  DO i = 1, nlines
    READ(1,*) energies(i,:)
  END DO
  CLOSE (1)

  DO i = 1, nlines
    energies(i,1) = energies(i,1) - deltaE
  END DO

  ! ..........................................................................!
  ! Initial calculation to establish baseline variance
  DO i = 1, nlines
    IF (energies(i,3) .NE. 0.0) THEN
      sum1 = sum1 + (energies(i,1)-energies(i,3))**2
      sum2 = sum2 + (energies(i,1)-energies(i,3))
    END IF
  END DO

  PRINT *, 'Delta E calc: ', (1/float(nlines)) * sum2

  sum1  = (1/float(nlines)) * sum1
  sum2  = (1/float(nlines)**2) * sum2**2
  total = sum1 - sum2
  criteria = float(n)*sqrt(total)

  PRINT *, 'Variance:      ', sqrt(total)

  ! End calculation of baseline variance
  ! ..........................................................................!

  ! ..........................................................................!
  ! Calculation for discarding data based on calculated variance
  PRINT *, 'Criteria: ', criteria

  DO i = 1, nlines
    compare = ABS(ABS(energies(i,1))-ABS(energies(i,3)))
    IF (compare .GT. criteria) THEN
      energies(i,:) = 0.0
      discard = discard + 1
    END IF
  END DO

  PRINT *, '# discarded data points: ', discard

  ! Add code to print shifted energies to file for Python plot
  ! Will need to modify when we start removing data
  OPEN(2,file = outFile)
  DO i = 1, nlines
    WRITE(2,*) energies(i,:)
  END DO
  CLOSE(2)

  ! Add ability to save multiple iterations of data to file after
  ! adding the minimization

END PROGRAM rmsFortran
