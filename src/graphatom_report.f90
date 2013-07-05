MODULE Graphatom_report
  USE atomdata
  USE gridmod
  USE report_mod

  IMPLICIT NONE

CONTAINS
  SUBROUTINE Report_Graphatom(key,Grid,Orbit,Pot,SCF)

    CHARACTER(2) :: key
    TYPE(GridInfo) :: Grid
    TYPE(OrbitInfo) :: Orbit
    TYPE(PotentialInfo) :: Pot
    TYPE(SCFInfo) :: SCF

    INTEGER :: i,j,io,many,l,istart
    CHARACTER (len=4) :: flnm
    CHARACTER (len=20) :: nm
    CHARACTER (len=2) :: sym
    CHARACTER (len=1) :: syml
    REAL(8), allocatable :: wav(:,:)
    REAL(8), POINTER :: r(:),den(:),rv(:),wfn(:,:)
    INTEGER :: n,nps,npp,npd,npf,npg
    INTEGER, PARAMETER :: ifen=10, ifden=11,ifwfn=12

    OPEN (unit = ifen, file=TRIM(Pot%sym)//'.GA', form='formatted', &
&        position = 'append')

    call summary_report(ifen,key,Grid,Orbit,Pot,SCF)

    CLOSE(ifen)

    write(6,*) 'finished summary'; call flush(6)
    IF (key=='AE') THEN
       n=Grid%n
       r=>Grid%r
       den=>Orbit%den
       rv=>Pot%rv
       nps=Orbit%nps
       npp=Orbit%npp
       npd=Orbit%npd
       npf=Orbit%npf
       npg=Orbit%npg
       wfn=>Orbit%wfn
       allocate(wav(n,Orbit%norbit))
       wav=wfn

       !
       ! write density and wavefunctions
       !

       OPEN (unit = ifden, file = 'density.GA', form = 'formatted')
       DO i = 1, n
          IF (r (i) .LE.6.d0) THEN
             WRITE (ifden, '(1p,6e12.4)') r (i), den (i)
          ENDIF
       ENDDO
       CLOSE (ifden)

       OPEN (unit = ifden, file = 'potential.GA', form = 'formatted')
       DO i = 1, n
          IF (r (i) .LE.6.d0) THEN
             WRITE (ifden, '(1p,6e12.4)') r (i), rv (i)
          ENDIF
       ENDDO
       CLOSE (ifden)

       OPEN (unit = ifden, file = 'plotdensity.GA', form = 'formatted')
       nm = 'density.GA'
       WRITE (ifden, '("gplot -t ""Radial density for ",a2, &
&             """ -tx ""r (bohr)"" -f ",a10, " 1 2 lines")') pot%sym, nm
       CLOSE (ifden)
       OPEN (unit = ifden, file = 'plotpotential.GA', form = 'formatted')
       nm = 'potential.GA'
       WRITE (ifden, '("gplot -t ""rxV(r) for ",a2, &
&             """ -tx ""r (bohr)"" -f ",a12, " 1 2 lines")') pot%sym, nm
       CLOSE (ifden)
       DO i = 1, n
          DO io = 1,Orbit%norbit
             IF (dabs (wav (i, io) ) .LT.1.d-8) wav (i, io) = 0.d0
          ENDDO
       ENDDO
       istart = 0
       DO l = 0, 4
          IF (l.EQ.0) many = nps
          IF (l.EQ.1) many = npp - 1
          IF (l.EQ.2) many = npd-2
          IF (l.EQ.3) many = npf - 3
          IF (l.EQ.4) many = npg - 4
          IF (l.EQ.0) syml = 's'
          IF (l.EQ.1) syml = 'p'
          IF (l.EQ.2) syml = 'd'
          IF (l.EQ.3) syml = 'f'
          IF (l.EQ.4) syml = 'g'
          IF (many.GT.0) THEN
             CALL mkname (l, flnm)
             OPEN (unit = ifwfn, file = 'GAwfn'//flnm, form = 'formatted')
             DO i = 1, n
                IF (r (i) .LE.6.d0) THEN
                   WRITE (ifwfn, '(1p,8e10.2)') r (i) , (wav (i, j) , j = &
&                       istart + 1, istart + many)
                ENDIF
             ENDDO
             CLOSE (ifwfn)
             nm = 'GAwfn'//flnm
             OPEN (unit = ifwfn, file = 'plotGAwfn'//flnm, form = 'formatted')
             IF (many.EQ.1) WRITE (ifwfn, '("gplot -t ""Radial ",a1, &
                  &   " wavefunctions for ",a2, &
                  &""" -tx ""r (bohr)"" -f ",a8, " 1 2 lines")') syml, pot%sym, nm
             IF (many.EQ.2) WRITE (ifwfn, '("gplot -t ""Radial ",a1, &
                  &   " wavefunctions for ",a2, &
                  &""" -tx ""r (bohr)"" -f ",a8, " 1 2 lines -f ",a8," 1 3 lines")') syml &
                  , pot%sym, nm, nm
             IF (many.EQ.3) WRITE (ifwfn, '("gplot -t ""Radial ",a1, &
                  &   " wavefunctions for ",a2, &
                  &""" -tx ""r (bohr)"" -f ",a8, " 1 2 lines -f ",a8," 1 3 lines -f ", &
                  &a8," 1 4 lines")') syml, pot%sym, nm, nm, nm
             IF (many.EQ.4) WRITE (ifwfn, '("gplot -t ""Radial ",a1, &
                  &   " wavefunctions for ",a2, &
                  &""" -tx ""r (bohr)"" -f ",a8, " 1 2 lines -f ",a8," 1 3 lines -f ", &
                  &a8," 1 4 lines -f ",a8," 1 5 lines")') syml, pot%sym, nm, nm, nm, nm
             IF (many.EQ.5) WRITE (ifwfn, '("gplot -t ""Radial ",a1, &
                  &   " wavefunctions for ",a2, &
                  &""" -tx ""r (bohr)"" -f ",a8, " 1 2 lines -f ",a8," 1 3 lines -f ", &
                  &a8," 1 4 lines -f ",a8," 1 5 lines -f ",a8," 1 6 lines")') syml,  &
                  &pot%sym, nm, nm, nm, nm, nm
             IF (many.EQ.6) WRITE (ifwfn, '("gplot -t ""Radial ",a1, &
                  &   " wavefunctions for ",a2, &
                  &""" -tx ""r (bohr)"" -f ",a8, " 1 2 lines -f ",a8," 1 3 lines -f ", &
                  &a8," 1 4 lines -f ",a8," 1 5 lines -f ",a8," 1 6 lines -f " &
                  &,a8," 1 7 lines")') syml, pot%sym, nm, nm, nm, nm, nm, nm
             IF (many.EQ.7) WRITE (ifwfn, '("gplot -t ""Radial ",a1, &
                  &   " wavefunctions for ",a2, &
                  &""" -tx ""r (bohr)"" -f ",a8, " 1 2 lines -f ",a8," 1 3 lines -f ", &
                  &a8," 1 4 lines -f ",a8," 1 5 lines -f ",a8," 1 6 lines -f " &
                  &,a8," 1 7 lines -f ",a8," 1 8 lines" )') syml, pot%sym, nm, nm, nm, nm &
                  &, nm, nm, nm
             CLOSE (ifwfn)
          ENDIF
          istart = istart + many

       ENDDO
    deallocate(wav)
    ENDIF
  END SUBROUTINE Report_Graphatom

END MODULE Graphatom_report
