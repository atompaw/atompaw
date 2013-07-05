!******************************************************************************
!
! File : search_sort.f90
!   by : Alan Tackett
!   on : 11/17/95
!  for : General Routines for Searching and Sorting arrays
!
! Searching Routines
!      Linear_Search - Performs a linear search
!      Binary_Search - Performs a binary search
!
! Sorting Routines
!      Insertion_Sort - Sorts an array using the Insertion Sort method
!
!
!******************************************************************************

MODULE search_sort

  IMPLICIT NONE

  INTERFACE Insertion_Sort
     MODULE PROCEDURE Real_InsSort
     MODULE PROCEDURE Integer_InsSort
  END INTERFACE


CONTAINS

  !******************************************************************************
  !
  !  Linear_Search - Performs a linear search of an array and returns
  !                  the index of the closest matching element
  !
  !  A          - Array of values to search
  !  MatchValue - Value to match
  !  Exact      - Upon return it contains TRUE if an exact match was found
  !               and FALSE otherwise. (OPTIONAL)
  !
  !  Return Values
  !       The closest match, index is returned. If EXACT is
  !       supplied then it determines whether an exact match was found.
  !
  !******************************************************************************

  INTEGER FUNCTION Linear_Search(A, MatchValue, Exact)
    INTEGER,           INTENT(IN)  :: A(:)
    INTEGER,           INTENT(IN)  :: MatchValue
    LOGICAL, OPTIONAL, INTENT(OUT) :: Exact

    LOGICAL  :: Found
    INTEGER  :: i, BestIndex, A_Size
    INTEGER  :: BestValue, BestMiss

    A_Size = SIZE(A)
    BestValue = A(1)
    BestIndex = 1
    BestMiss = ABS(BestValue - MatchValue)
    i = 2
    Found = .FALSE.

    DO WHILE ((i<=A_Size) .AND. (BestValue /= MatchValue))
       IF (ABS(A(i) - MatchValue) < BestMiss) THEN
          BestIndex = i
          BestValue = A(i)
          BestMiss = ABS(BestValue - MatchValue)
       END IF

       i = i + 1
    END DO

    IF (PRESENT(Exact)) THEN
       IF (BestValue == MatchValue) THEN
          Exact = .TRUE.
       ELSE
          Exact = .FALSE.
       END IF
    END IF

    Linear_Search = BestIndex
    RETURN
  END FUNCTION Linear_Search


  !******************************************************************************
  !
  !  Linear_Search - Performs a binary search of an array and returns
  !                  the index of the closest matching element that
  !                  is LESS THAN OR EQUAL to the element.
  !
  !  A          - Array of values to search
  !  MatchValue - Value to match
  !  Exact      - Upon return it contains TRUE if an exact match was found
  !               and FALSE otherwise. (OPTIONAL)
  !
  !  Return Values
  !       The closest match, index is returned. If EXACT is
  !       supplied then it determines whether an exact match was found.
  !
  !******************************************************************************

  INTEGER FUNCTION Binary_Search(A, MatchValue, Exact)
    INTEGER,           INTENT(IN)  :: A(:)
    INTEGER,           INTENT(IN)  :: MatchValue
    LOGICAL, OPTIONAL, INTENT(OUT) :: Exact

    INTEGER  :: A_Size
    INTEGER  :: HiIndex, LoIndex, i

    A_Size = SIZE(A)
    IF (A(1) > A(A_Size)) THEN
       HiIndex = 1
       LoIndex = A_Size
    ELSE
       HiIndex = A_Size
       LoIndex = 1
    END IF

    i = (HiIndex + LoIndex) / 2
    IF (A(HiIndex) <= MatchValue) THEN
       i = HiIndex
       LoIndex = i
    ELSE IF (A(LoIndex) >= MAtchValue) THEN
       i = LoIndex
       HiIndex = i
    END IF

    DO WHILE ((HiIndex /= LoIndex) .AND. (A(i) /= MAtchValue))
       IF (A(i) > MatchValue) THEN
          HiIndex = i
       ELSE
          LoIndex = i
       END IF

       i = (HiIndex + LoIndex) / 2
    END DO

    IF (PRESENT(Exact)) Exact = .FALSE.

    IF (A(i) == MatchValue) THEN
       IF (PRESENT(Exact)) Exact = .TRUE.
    END IF

    Binary_Search = i
    RETURN
  END FUNCTION Binary_Search

  !******************************************************************************
  !
  !  Insertion_Sort - Sorts an Array using an Insertion sort method.
  !                   Best used for small arrrays or arrays that are almost
  !                   sorted.
  !
  !  A         - Array to Sort
  !  LUT       - Index table containing the new sorted order(RETURNED)
  !  Ascending - Sort Direction.  If TRUE then sorted in ascending order
  !              otherwise the array is sorted in descending order.
  !
  !
  !******************************************************************************

  SUBROUTINE Integer_InsSort(A, LUT, Ascending)
    INTEGER, INTENT(IN)    :: A(:)
    INTEGER, INTENT(INOUT) :: LUT(:)
    LOGICAL, INTENT(IN)    :: Ascending

    INTEGER :: i, j, k, A_Size
    INTEGER :: CurrVal

    A_Size = SIZE(A)

    DO i=1, A_Size
       LUT(i) = i
    END DO

    DO i = 2, A_Size
       CurrVal = A(LUT(i))
       j = i - 1

       DO WHILE ((CurrVal < A(LUT(j))) .AND. (j>1))
          LUT(j+1) = LUT(j)
          j = j - 1
       END DO

       IF (CurrVal < A(LUT(j))) THEN
          LUT(j+1) = LUT(j)
          j = j - 1
       END IF

       LUT(j+1) = i
    END DO

    IF (.NOT. Ascending) THEN
       j = A_Size / 2

       DO i = 1, j
          k = LUT(i)
          LUT(i) = LUT(A_Size - i + 1)
          LUT(A_Size - i + 1) = k
       END DO
    END IF

    RETURN
  END SUBROUTINE Integer_InsSort

  !*******

  SUBROUTINE Real_InsSort(A, LUT, Ascending)
    REAL(8),    INTENT(IN)    :: A(:)
    INTEGER, INTENT(INOUT) :: LUT(:)
    LOGICAL, INTENT(IN)    :: Ascending

    INTEGER :: i, j, k, A_Size
    REAL(8)    :: CurrVal

    A_Size = SIZE(A)

    DO i=1, A_Size
       LUT(i) = i
    END DO

    DO i = 2, A_Size
       CurrVal = A(LUT(i))
       j = i - 1

       DO WHILE ((CurrVal < A(LUT(j))) .AND. (j>1))
          LUT(j+1) = LUT(j)
          j = j - 1
       END DO

       IF (CurrVal < A(LUT(j))) THEN
          LUT(j+1) = LUT(j)
          j = j - 1
       END IF

       LUT(j+1) = i
    END DO

    IF (.NOT. Ascending) THEN
       j = A_Size / 2

       DO i = 1, j
          k = LUT(i)
          LUT(i) = LUT(A_Size - i + 1)
          LUT(A_Size - i + 1) = k
       END DO
    END IF

    RETURN
  END SUBROUTINE Real_InsSort

END MODULE search_sort


