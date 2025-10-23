!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  This module is used only by the program aeatom and has 
!    the single subroutine:  Report_aeonly
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

MODULE aeonly_report

  USE Tools
  USE atomdata
  USE gridmod
  USE report_mod

  IMPLICIT NONE

CONTAINS
  SUBROUTINE Report_aeonly(key,Grid,Orbit,Pot,SCF)

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
    REAL(8), allocatable :: wav(:,:),lwav(:,:)
    REAL(8), POINTER :: r(:),den(:),rv(:),wfn(:,:),lwfn(:,:)
    INTEGER :: n,nps,npp,npd,npf,npg
    INTEGER, PARAMETER :: ifen=10, ifden=11,ifwfn=12

    OPEN (unit = ifen, file=TRIM(Pot%sym)//'.AEO', form='formatted', &
&        position = 'append')

    call summary_report(ifen,key,Grid,Orbit,Pot,SCF)

    CLOSE(ifen)

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
       lwfn=>Orbit%lwfn

       allocate(wav(n,Orbit%norbit))
       if (diracrelativistic) allocate(lwav(n,Orbit%norbit))
       wav=wfn
       if (diracrelativistic) lwav=lwfn

       !
       ! write density and wavefunctions
       !

       OPEN (unit = ifden, file = 'density.AEO', form = 'formatted')
       DO i = 1, n
          IF (r (i) .LE.6.d0) THEN
             WRITE (ifden, '(1p,6e12.4)') r (i), den (i)
          ENDIF
       ENDDO
       CLOSE (ifden)

       OPEN (unit = ifden, file = 'potential.AEO', form = 'formatted')
       DO i = 1, n
          IF (r (i) .LE.6.d0) THEN
             WRITE (ifden, '(1p,6e12.4)') r (i), rv (i)
          ENDIF
       ENDDO
       CLOSE (ifden)

       OPEN (unit = ifden, file = 'plotdensity.AEO', form = 'formatted')
       nm = 'density.AEO'
       WRITE (ifden, '("gplot -t ""Radial density for ",a2, &
&             """ -tx ""r (bohr)"" -f ",a10, " 1 2 lines")') pot%sym, nm
       CLOSE (ifden)
       OPEN (unit = ifden, file = 'plotpotential.AEO', form = 'formatted')
       nm = 'potential.AEO'
       WRITE (ifden, '("gplot -t ""rxV(r) for ",a2, &
&             """ -tx ""r (bohr)"" -f ",a12, " 1 2 lines")') pot%sym, nm
       CLOSE (ifden)
       DO i = 1, n
          DO io = 1,Orbit%norbit
             IF (dabs (wav (i, io) ) .LT.1.d-8) wav (i, io) = 0.d0
             if (diracrelativistic) then 
              IF (dabs (lwav (i, io) ) .LT.1.d-8) lwav (i, io) = 0.d0
             endif 
          ENDDO
       ENDDO
       if (.not.diracrelativistic) then
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
             OPEN (unit = ifwfn, file = 'AEOwfn'//flnm, form = 'formatted')
             DO i = 1, n
                IF (r (i) .LE.6.d0) THEN
                   WRITE (ifwfn, '(1p,100e10.2)') r (i) , (wav (i, j) , j = &
&                       istart + 1, istart + many)
                ENDIF
             ENDDO
             
             CLOSE (ifwfn)
          ENDIF
          istart = istart + many

       ENDDO
       endif
       if(diracrelativistic) then
          OPEN(unit=ifwfn, file='AEOwfnall', form='formatted')     
             DO i = 1, n
                IF (r (i) .LE.6.d0) THEN
                   WRITE (ifwfn, '(1p,100e10.2)') r (i) , (wav (i, j) ,&
&                     lwav(i,j),  j = 1,Orbit%norbit)
                ENDIF
             ENDDO
          CLOSE(ifwfn)     
       endif   
    deallocate(wav)
    if (diracrelativistic) deallocate(lwav)
    ENDIF

  END SUBROUTINE Report_aeonly

END MODULE aeonly_report
