!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  This module contains the following active subroutines:
!      mkname, PrintDate, PrintDateStr, ConvertChar, flush_unit
!      extractword, UpperCase, eliminate_comment
!  This module contains the following active functions:
!      stripchar, checkline2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

MODULE Tools

  IMPLICIT NONE

CONTAINS

!******************************************************************************
!  mkname - Subroutine to take an integer (.le. 4 digits) and return it
!           in the form of a character string
!******************************************************************************
  SUBROUTINE mkname(i,stuff)
    CHARACTER(4) stuff
    INTEGER i,i1,i10,i100,i1000
    stuff='?'
    IF (i.GT.9999) i=MOD(i,10000)
    i1000=i/1000
    i100=(i-1000*i1000)/100
    i10=(i-1000*i1000-100*i100)/10
    i1=(i-1000*i1000-100*i100-10*i10)
    IF (i.GE.1000) THEN
       stuff=CHAR(i1000+48)//CHAR(i100+48)//CHAR(i10+48)//CHAR(i1+48)
       RETURN
    ENDIF
    IF (i.GE.100) THEN
       stuff=CHAR(i100+48)//CHAR(i10+48)//CHAR(i1+48)
       RETURN
    ENDIF
    IF (i.GE.10) THEN
       stuff=CHAR(i10+48)//CHAR(i1+48)
       RETURN
    ENDIF
    IF (i.GE.0) stuff=CHAR(i1+48)
    RETURN
  END SUBROUTINE mkname

!******************************************************************************
!  PrintDate - Prints the date to the specified unit with TEXT prepended.
!    Unit - Output unit
!    Text - Text for prepending. All trailing blanks are removed.
!******************************************************************************
  SUBROUTINE PrintDate(Unit, Text)
    INTEGER,       INTENT(IN) :: Unit
    CHARACTER*(*), INTENT(IN) :: Text
    CHARACTER*10 :: DateStr, TimeStr
    CHARACTER*50 :: FmtStr
    CALL DATE_AND_TIME(DateStr, TimeStr)
    FmtStr=' ' // DateStr(5:6) // '/' // DateStr(7:8) // '/' // DateStr(1:4) // &
&        ', ' // TimeStr(1:2) // ':' // TimeStr(3:4) // ':' // &
&        TimeStr(5:10)
    WRITE(Unit, '(a,a)')  TRIM(Text), TRIM(FmtStr)
    RETURN
  END SUBROUTINE PrintDate

!******************************************************************************
!  PrintDateStr - Prints the date to the specified string
!******************************************************************************
Subroutine PrintDateStr(text)
  Character*(*), Intent(out) :: Text
  Character*10 :: DateStr, TimeStr
  Character*50 :: FmtStr
  Call Date_And_Time(DateStr, TimeStr)

  FmtStr=' ' // DateStr(5:6) // '/' // DateStr(7:8) // '/' // DateStr(1:4) // &
       &       ', ' // TimeStr(1:2) // ':' // TimeStr(3:4) // ':' // &
       &         TimeStr(5:10)
  Write(text, '(a)')  Trim(FmtStr)
  Return
End Subroutine

!******************************************************************************
!  stripchar - Eliminate blanks
!******************************************************************************
  FUNCTION stripchar(inputchar)
    CHARACTER(132) :: stripchar
    CHARACTER*(*), INTENT(IN) :: inputchar
    INTEGER :: i,j,n
    n=LEN(inputchar)
    DO i=1,132
       stripchar(i:i)=''
    ENDDO
    j=0
    DO i=1,n
       IF (inputchar(i:i) /= '') THEN
          j=j+1
          stripchar(j:j)=inputchar(i:i)
       ENDIF
    ENDDO
  END FUNCTION stripchar

!******************************************************************************
!  ConvertChar - convert a string into an integer
!******************************************************************************
  SUBROUTINE ConvertChar(inchar,outn)
    CHARACTER*(*),INTENT(IN) :: inchar
    INTEGER, INTENT(OUT) :: outn
    INTEGER :: i, j, k, n, fac
    n=LEN(inchar)
    fac=1; outn=0
    DO i=n,1,-1
       IF (inchar(i:i)==''.OR.inchar(i:i)=="-") EXIT
       j=ICHAR(inchar(i:i))-48
       outn=outn+fac*j
       fac=fac*10
    ENDDO
    outn=outn
  END SUBROUTINE ConvertChar

!******************************************************************************
!  extractword - extract a word from a string
!******************************************************************************
  SUBROUTINE extractword(wordindex,stringin,stringout)
    INTEGER, INTENT(IN) :: wordindex
    CHARACTER(*), INTENT(IN) :: stringin
    CHARACTER(*), INTENT(OUT) :: stringout
    INTEGER :: i,j,n,str,fin,icount
    stringout=''
    n=LEN(stringin)
    i=INDEX(stringin,'!');IF (i==0) i=n
    j=INDEX(stringin,'#');IF (j==0) j=n
    n=MIN(i,j,n)
    str=1;fin=n
    DO icount=1,MAX(1,wordindex-1)
       DO i=str,n
          IF (stringin(i:i)/=' ') EXIT
       ENDDO
       str=i
       IF (n>str) THEN
          DO i=str+1,n
             IF(stringin(i:i)==' ') EXIT
          ENDDO
          fin=i
       ENDIF
       IF (wordindex>2) THEN
          IF (fin<n) THEN
             str=fin+1
          ELSE
             EXIT
          ENDIF
       ENDIF
    ENDDO
    IF (wordindex>1) THEN
       IF (fin>=n) RETURN
       DO i=fin+1,n
          IF (stringin(i:i)/=' ') EXIT
       ENDDO
       str=i
       IF (n>str) THEN
          DO i=str+1,n
             IF(stringin(i:i)==' ') EXIT
          ENDDO
          fin=i
       ENDIF
    ENDIF
    stringout=stringin(str:fin)
  END SUBROUTINE extractword

!******************************************************************************
!  UpperCase - Converts a string to Upper Case
!  str     - String to convert
!******************************************************************************
  SUBROUTINE UpperCase(str)
    CHARACTER(*), INTENT(INOUT) :: str
    INTEGER  :: i, j, k
    j = LEN(Str)
    DO i=1, j
       k = IACHAR(str(i:i))
       IF ((k>96) .AND. (k<123)) str(i:i) = ACHAR(k-32)
    END DO
    RETURN
  END SUBROUTINE UpperCase

!******************************************************************************
!  checkline2 - logical function
!     Returns .true. if inputline matchs either in1 or in2
!     Note: in1 and in2 should either be integers or characters
!******************************************************************************
  LOGICAL FUNCTION checkline2(inputline,in1,in2)
     CHARACTER(*), INTENT(IN) :: inputline
     CHARACTER(*), INTENT(IN) :: in1,in2
     INTEGER :: i,leninput,len1,len2
     CHARACTER(120) :: inputline_u,in1_u,in2_u
     inputline_u=trim(inputline) ; call UpperCase(inputline_u)
     in1_u=trim(in1) ; call UpperCase(in1_u)
     in2_u=trim(in2) ; call UpperCase(in2_u)
     leninput=len(trim(inputline));len1=len(trim(in1));len2=len(trim(in2))
     checkline2=.false.
     if (leninput==len1) checkline2=(inputline_u(1:len1)==trim(in1))
     if ((.not.checkline2).and.leninput==len2) checkline2=(inputline_u(1:len2)==trim(in2))
     RETURN
  END FUNCTION checkline2

!******************************************************************************
!  eliminate_comment - Eliminate comment on the right of a line (! or #)
!     line - string to convert (output replaces input)
!******************************************************************************
  SUBROUTINE eliminate_comment(line)
    CHARACTER(*), INTENT(INOUT) :: line
    INTEGER :: i,i0
    i0=-1 ; i=1
    DO WHILE (i0<0.AND.i<LEN(line))
      i=i+1
      IF (line(i:i)=="!".OR.line(i:i)=="#") i0=i
    END DO
    IF (i0 >1) line=line(1:i0-1)
    IF (i0==1) line=""
  END SUBROUTINE eliminate_comment

END MODULE Tools
