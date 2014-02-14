C
C ----------------------------------------------------------------
C     Subroutine RANECU
C
C     Portable random number generator proposed by l'Ecuyer
C                 in Commun. ACM 31 (1988) 743.
C     slightly modified by F. James, 1988, to generate a vector
C     of pseudorandom numbers RVEC of length LEN.
C     Listing Published in Computer Physics COmm. 60 (1990) 329.
C     Typed in by S. A. Wood, CEBAF 1/23/91 with strong typing.
C ----------------------------------------------------------------
C
      SUBROUTINE RANECU(RVEC,LEN)
      IMPLICIT NONE
C
      REAL RVEC(*)
      INTEGER LEN,ISEED1,ISEED2
      INTEGER I,K,IZ,IS1,IS2
C
      SAVE ISEED1,ISEED2
      DATA ISEED1,ISEED2/12345,67890/
C
      DO I=1,LEN
        K = ISEED1/53668
        ISEED1 = 40014*(ISEED1-K*53668)-K*12211
        IF(ISEED1.LT.0) ISEED1 = ISEED1+2147483563
C
        K = ISEED2/52774
        ISEED2 = 40692*(ISEED2-K*52774)-K*3791
        IF(ISEED2.LT.0) ISEED2 = ISEED2+2147483399
C
        IZ = ISEED1 - ISEED2
        IF(IZ.LT.1) IZ = IZ + 2147483562
C
        RVEC(I) = REAL(IZ)*4.656613E-10
      ENDDO
      RETURN
C
C     Initialize seeds.
C
      ENTRY RECUIN(IS1,IS2)
      ISEED1 = IS1
      ISEED2 = IS2
      RETURN
C
C     Get current seeds.
C
      ENTRY RECUUT(IS1,IS2)
      IS1 = ISEED1
      IS2 = ISEED2
      RETURN
C
      END
