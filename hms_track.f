*------------------------------------------------------------------------
*
*       HMS_TRACK  HMS Tracking Routines 
*      -=========-
* 
*	Forward and Backward Tracking of electrons in the Jlab HMS hall 
*       C spectrometer
*
*       Note: - the HMS routines use a lab (HMS) coord. system
*               and the corresponding COSY coord. system, both 
*               right handed with
*                 x : pointing downwards
*                 y : perpendicular to x,z, 
*                     pointing to the left (if seen in z-direction)
*                 z : HMS axis, pointing from the target to the focal plane
*
*             - all lengths (x,y,z,l,...) are measured in [m]
*             - all angles are measured as dx/dz,dy/dz (lab coords.)
*               or as A,B (COSY coords.)       
*             - the momentum is measured in delta (relative momentum
*               deviation  = 1-p0/pHMS)
*
*       PART 1: Forward trackinh using COSY transport matrices  
*              
*
*       PART 2: Reconstruction (backward tracking) using reconstruction 
*               and COSY transport matrices (including the 
*               effects of a vertical beam offset (out-of plane))  
*                    
*
*       written by Markus Muehlbauer for the GEN Experiment
*
* frw 9/2000
*       changes made in the course of the migration to g77 compiler:
*        - fixed some typos (old ones, too)
*        - all COSY conversion code was already commented out (why?
*          by who?) so I removed it and made the code more readable
*        - various variables were not initialized prior to reading
*          from file.  This may or may not be an issue, but fixed it
*          anyway
*
*------------------------------------------------------------------------

*------------------------------------------------------------------------
*
*       PART 1:  HMS Forward Tracking (Target to Focal Plane)
*      -=======-       
*
*	Forward tracking in the Jlab HMS hall C spectrometer
*       using COSY transport matrices
*   
*       developed by Cris Cothran
*       modified  by Markus Muehlbauer
*         - CCs orignal program converted into subroutines
*         - mad additions for pure tracking, without checking the acceptance
*         - and changed the innermost loops applying the matrix
*           (which speeds up the whole thing by a factor of about 30)
* 
*       Supplies:
*         hmsInitForward (map) 
*           load the forward transport maps
*         hmsForward (uT,zT,u,z)
*            make a single step transport calculation 
*            (without treating the acceptance)
*         hmsAccept (uT,zT,u,z)
*            make a multi step transport calculation 
*            (also treating the acceptance)  
*       
*       Note: - Before calling hmsForward or hmsAccept the forward 
*               transport maps have to be loaded by a call to hmsInitForward
*------------------------------------------------------------------------

      SUBROUTINE hmsInitForward (filename,maxorder,path,p0) 
      IMPLICIT  NONE
      CHARACTER filename*(*)
      INTEGER   maxorder
      LOGICAL   path
      REAL*8      p0
      
* --  load the HMS forward (cosy) map 
*
*     Parameter:
*       filename   I : name of the forward map
*       maxorder   I : maximal order to take into account
*       path       I : calculate the path length/TOF variation
*       p0         I : momentum the spectrometer is set to
*
*     reconInitFormat uses the same map-file as the simulation
*     routines and digs out the necessary matrix elements to
*     perform the focal plane corrections
* 
*     Map Line Ax Ax' Ay Ay' Al Adelta XxYyld
*     
*       Map    : Map number 0-20
*                20: full HMS (target to focal plane)
*       Line   : Line number
*       Ax     : coeff's applied to x 
*       Ax'    : coeff's applied to x'
*       Ay     : coeff's applied to y
*       Ay'    : coeff's applied to y'
*       Al     : coeff's applied to path length deviation
*       X      : exponent of the actual x
*       x      : exponent of the actual x'
*       Y      : exponent of the actual y
*       y      : exponent of the actual y'
*       l      : exponent of the actual path length
*       d      : exponent of the delta

      REAL*8       fac,me
      PARAMETER (fac = 0.9904599d00)
      PARAMETER (me  = 0.00051099906d00)

      include 'trans_map.inc'

      ! other variables
      INTEGER   i,m,n,num,eof 
       
      ! read the necessary coeffs from the map data file
 10   FORMAT (1X,I2,1X,I5,1X,5(E14.7E2,1X),6I1)
 
      mep0 = me/(p0*fac)
 
      DO i=1,NMAPS
        first(i) = 0
        last (i) = 0 
      ENDDO
 
      DO i=1,NLINES
       c1(i) = 0.d00 
       c2(i) = 0.d00 
       c3(i) = 0.d00 
       c4(i) = 0.d00 
       c5(i) = 0.d00        
       e1(i) = 0 
       e2(i) = 0 
       e3(i) = 0 
       e4(i) = 0 
       e5(i) = 0 
       e6(i) = 0 
      ENDDO

      num   = 1
      order = 0

      OPEN (UNIT=97,FILE=filename,STATUS='OLD')
      
      READ (97,10,IOSTAT=eof)  m,n, 
     >  c1(num), c2(num), c3(num), c4(num), c5(num), 
     >  e1(num), e2(num), e3(num), e4(num), e5(num), e6(num) 
      IF (.not. path) c5(num) = 0.
       
      DO WHILE (eof .GE. 0)
        
        i = e1(num)+e2(num)+e3(num)+e4(num)+e5(num)+e6(num)
        IF ((i .LE. maxorder) .AND. ((c1(num) .NE. 0.0d00) .OR.
     >       (c2(num) .NE. 0.0) .OR. (c3(num) .NE. 0.0d00) .OR.
     >       (c4(num) .NE. 0.0) .OR. (c5(num) .NE. 0.0d00))) THEN
          IF (i .GT. order)    order    = i
          IF (0 .EQ. first(m)) first(m) = num
          last(m) = num
          num     = num+1
        ENDIF
        READ (97,10,IOSTAT=eof)  m,n, 
     >    c1(num), c2(num), c3(num), c4(num), c5(num), 
     >    e1(num), e2(num), e3(num), e4(num), e5(num), e6(num) 
        IF (.not. path) c5(num) = 0.
        
      ENDDO
      
      CLOSE (97)
 
      RETURN
      END 



*------------------------------------------------------------------------

      SUBROUTINE hmsApplyCOSY (u,num)
      IMPLICIT NONE
      REAL*8    u(6)
      INTEGER num
      
* --  apply a COSY matrix on the COSY vector u 
*     (used in the forward tracking only)
*
*     Parameter
*       u   IO : coordinate vector (COSY)
*                 u(1,2) : x [m], A = out of plane coordinates (downwards) 
*                 u(3,4) : y [m], B = inplane coordinates (perp. on x,z)
*                 u(5)   : l [m] = [?] = path length deviation (?)
*                 u(6)   : delta = relative deviation of the particle 
*                                  momentum from p0
*       num I  : cosy matrix to apply 

      include 'trans_map.inc'

      ! other variables
      REAL*8    a,uu1(0:10),uu2(0:10),uu3(0:10) 
      REAL*8      uu4(0:10),uu5(0:10),uu6(0:10)
      INTEGER i 
         
      ! calculate the powers of the focal plane coordinates
      uu1(0) = 1.d00
      uu2(0) = 1.d00
      uu3(0) = 1.d00
      uu4(0) = 1.d00
      uu5(0) = 1.d00
      uu6(0) = 1.d00
      uu1(1) = u(1)
      uu2(1) = u(2)
      uu3(1) = u(3)
      uu4(1) = u(4)
      uu5(1) = u(5)
      uu6(1) = u(6)
      DO i=2,order
        uu1(i)=uu1(i-1)*uu1(1)
        uu2(i)=uu2(i-1)*uu2(1)
        uu3(i)=uu3(i-1)*uu3(1)
        uu4(i)=uu4(i-1)*uu4(1)
        uu5(i)=uu5(i-1)*uu5(1)
        uu6(i)=uu6(i-1)*uu6(1)
      ENDDO
     
      DO i=1,5
        u(i) = 0.d00
      ENDDO
      
      ! apply the cosy matrix
      DO i = first(num),last(num) 
        a = uu1(e1(i)) * uu2(e2(i)) * uu3(e3(i)) * 
     >      uu4(e4(i)) * uu5(e5(i)) * uu6(e6(i))
        u(1) = u(1) + c1(i) * a
        u(2) = u(2) + c2(i) * a
        u(3) = u(3) + c3(i) * a
        u(4) = u(4) + c4(i) * a
        u(5) = u(5) + c5(i) * a
      ENDDO
    
      RETURN
      END   

*------------------------------------------------------------------------

      SUBROUTINE hmsForward (uT,zT,u,zz,lost)
      IMPLICIT none
      REAL*8     uT(6),zT,u(6),zz
      LOGICAL  lost
      
* --  make a single step transport calculation (without treating the acceptance)
*
*     Parameter:
*       uT,zT  I : target coordinates
*                    uT(1,2) : x [m], dx/dz = out of plane coord. (downwards) 
*                    uT(3,4) : y [m], dy/dz = inplane coord. (perp. on x,z)
*                    uT(5)   : l [?] = path length deviation (?) 
*                    uT(6)   : delta; relative deviation of the particle 
*                                     momentum from p0 
*                    zT      : z [m]; in axis coordinate (towards HMS)
*       u,zz   O : focal plane coordinates  
*                    u(1,2)  : x [m], dx/dz; out of plane coord. (downwards) 
*                    u(3,4)  : y [m], dy/dz; inplane coord. (perp. on x,z)
*                    u(5)    : l [?]; path length deviation (?) 
*                    u(6)    : delta; relative deviation of the particle 
*                                     momentum from p0 
*                    zz      : z [m]; in axis coordinate (towards HMS)
*                               position where the particle stops inside HMS
*                               (here always z-calorimeter) 
*       lost   O : set to .true. if the particle is lost in the HMS, 
*                  otherwise .false. (here always false)         

      ! other variables
      INTEGER i
    
      lost = .FALSE.
      zz   = 26.44743d00

      ! copy the target coordinates to the u vector 
      DO i=1,6
        u(i)=uT(i)
      ENDDO

      ! drift backwards to z=0 
      u(1)=u(1)-u(2)*zT
      u(3)=u(3)-u(4)*zT
 
      !go up to the hut
      CALL hmsApplyCOSY (u,20)
         
      RETURN
      END 

*------------------------------------------------------------------------

      LOGICAL FUNCTION hmsCheckDipole (u,zz,dz)
      IMPLICIT NONE
      REAL*8     u(6),zz,dz
      
* --  check for the dipole aperture
*
*     Parameter
*       u       I : target coordinates (COSY or LAB)
*                     u(1,2) : x [m], A or dx/dz = out of plane coord. (downwards) 
*                     u(3,4) : y [m], A or dy/dz = inplane coord. (perp. on x,z)
*                     u(5)   : l [?] = path length deviation (?) 
*                     u(6)   : delta; relative deviation of the particle 
*                                     momentum from p0 
*       zz     IO : z [m]; in axis coordinate (towards HMS)
*       dz      I : distance from last check point [m]
      
      zz = zz+dz
 
      hmsCheckDipole = .TRUE.
    
      IF (ABS(u(1)).GT.0.27940d00.OR.ABS(u(3)).GT.0.18415d00) THEN
        IF (ABS(u(1)).GT.0.34290d00.OR.ABS(u(3)).GT.0.20320d00) RETURN
        IF (ABS(u(1)).GT.0.27940d00.AND.ABS(u(3)).GT.0.12065d00) THEN
          IF (((ABS(u(1))-0.27940d00)**2 +
     >         (ABS(u(3))-0.12065d00)**2).GT.(0.06350d00)**2) RETURN
        ENDIF
        IF (ABS(u(1)).GT.0.13970d00 .OR.
     >     (ABS(u(1))-10.1852d00*ABS(u(3))).GT.2.069633d00) RETURN
      ENDIF
      
      hmsCheckDipole = .FALSE.
      
      RETURN
      END

*------------------------------------------------------------------------

      LOGICAL FUNCTION hmsCheckQuad (u,zz,dz,r)
      IMPLICIT NONE
      REAL*8     u(6),zz,dz,r
      
* --  check for the quadrupole aperture
*
*     Parameter
*       u       I : target coordinates (COSY or LAB)
*                     u(1,2) : x [m], A or dx/dz = out of plane coord. (downwards) 
*                     u(3,4) : y [m], A or dy/dz = inplane coord. (perp. on x,z)
*                     u(5)   : l [?] = path length deviation (?) 
*                     u(6)   : delta; relative deviation of the particle 
*                                     momentum from p0 
*       zz     IO : z [m]; in axis coordinate (towards HMS)
*       dz      I : distance from last check point [m]
*       r       I : aperture radius [m]
      
      zz = zz+dz
      hmsCheckQuad = ((u(1)**2+u(3)**2) .GT. r**2) 
           
      RETURN
      END

*------------------------------------------------------------------------

      LOGICAL FUNCTION hmsDriftOcta (u,zz,dz,x,y,m,b)
      IMPLICIT NONE
      REAL*8     u(6),zz,dz,x,y,m,b
      
* --  drift electron and check for the octagon
*
*     Parameter
*       u       I : target coordinates (COSY or LAB)
*                     u(1,2) : x [m], A or dx/dz = out of plane coord. (downwards) 
*                     u(3,4) : y [m], A or dy/dz = inplane coord. (perp. on x,z)
*                     u(5)   : l [?] = path length deviation (?) 
*                     u(6)   : delta; relative deviation of the particle 
*                                     momentum from p0 
*       zz     IO : z [m]; in axis coordinate (towards HMS)
*       dz      I : distance from last check point [m]
*       x,y,m,b I : octagon coordinates (x,y,b [m], m [1])

      zz=zz+dz
      u(1)=u(1)+u(2)*dz
      u(3)=u(3)+u(4)*dz
      hmsDriftOcta =  (ABS(u(1)) .GT. x) .OR. (ABS(u(3)) .GT. y) .OR.
     >               ((ABS(u(1))+m*ABS(u(3))) .GT. b)  
 
      RETURN
      END

*------------------------------------------------------------------------
      
      LOGICAL FUNCTION hmsDriftTPlate (u,zz,dz,r,y)
      IMPLICIT NONE
      REAL*8     u(6),zz,dz,r,y
      
* --  drift electron and check for the dipole transition plate
*
*     Parameter
*       u          IO : coordinate vector
*       zz         IO : z position [m]
*       dz         I  : drift distance
*       r,y        I  : aperture

      zz=zz+dz
      u(1)=u(1)+u(2)*dz
      u(3)=u(3)+u(4)*dz
      hmsDriftTPlate = ((u(1)**2+u(3)**2) .GT. r**2) .OR. 
     >                  (ABS(u(3)) .GT. y)
       
      RETURN
      END

*------------------------------------------------------------------------

      LOGICAL FUNCTION hmsDriftCirc (u,zz,dz,r)
      IMPLICIT NONE
      REAL*8     u(6),zz,dz,r

* --  drift electron and check for circular aperture
*
*     Parameter
*       u       I : target coordinates (COSY or LAB)
*                     u(1,2) : x [m], A or dx/dz = out of plane coord. (downwards) 
*                     u(3,4) : y [m], A or dy/dz = inplane coord. (perp. on x,z)
*                     u(5)   : l [?] = path length deviation (?) 
*                     u(6)   : delta; relative deviation of the particle 
*                                     momentum from p0 
*       zz     IO : z [m]; in axis coordinate (towards HMS)
*       dz      I : distance from last check point [m]
*       r,y     I : transition plate aperture radius [m]

      zz=zz+dz
      u(1)=u(1)+u(2)*dz
      u(3)=u(3)+u(4)*dz
      hmsDriftCirc = ((u(1)**2+u(3)**2) .GT. r**2)
       
      RETURN
      END

*------------------------------------------------------------------------
 
      LOGICAL FUNCTION hmsDriftRect (u,zz,dz,x0,x,y0,y)
      IMPLICIT NONE
      REAL*8     u(6),zz,dz,x0,x,y0,y 
      
* --  drift electron and check for rectangular aperture
*
*     Parameter
*       u       I : target coordinates (COSY or LAB)
*                     u(1,2) : x [m], A or dx/dz = out of plane coord. (downwards) 
*                     u(3,4) : y [m], A or dy/dz = inplane coord. (perp. on x,z)
*                     u(5)   : l [?] = path length deviation (?) 
*                     u(6)   : delta; relative deviation of the particle 
*                                     momentum from p0 
*       zz     IO : z [m]; in axis coordinate (towards HMS)
*       dz      I : distance from last check point [m]
*       x,x0,
*       y,y0    I : rectangular aperture [m] (half size and offset)

      zz=zz+dz
      u(1)=u(1)+u(2)*dz
      u(3)=u(3)+u(4)*dz
      hmsDriftRect = (ABS(u(1)-x0) .GT. x) .OR.
     >               (ABS(u(3)-y0) .GT. y) 
      RETURN
      END

*------------------------------------------------------------------------
            
      SUBROUTINE hmsAccept (uT,zT,u,zz,lost)
      IMPLICIT none
      REAL*8     uT(6),zT,u(6),zz
      LOGICAL  lost
      
* --  make a transport calculation to find the acceptance
*
*     Parameter:
*       uT,zT  I : target coordinates
*                    uT(1,2) : x [m], dx/dz; out of plane coord. (downwards) 
*                    uT(3,4) : y [m], dy/dz; inplane coord. (perp. on x,z)
*                    uT(5)   : l [?]; path length deviation (?) 
*                    uT(6)   : delta; relative deviation of the particle 
*                                     momentum from p0 
*                    zT      : z [m]; in axis coordinate (towards HMS)
*       u,zz   O : focal plane coordinates  
*                    u(1,2)  : x [m], dx/dz; out of plane coord. (downwards) 
*                    u(3,4)  : y [m], dy/dz; inplane coord. (perp. on x,z)
*                    u(5)    : l [?]; path length deviation (?) 
*                    u(6)    : delta; relative deviation of the particle 
*                                     momentum from p0 
*                    zz      : z [m]; in axis coordinate (towards HMS)
*                               position where the particle stops inside HMS
*       lost   O : set to .true. if the particle is lost in the HMS, 
*                  otherwise .false.          
 
      LOGICAL hmsCheckQuad 
      LOGICAL hmsCheckDipole 
      LOGICAL hmsDriftOcta
      LOGICAL hmsDriftTPlate
      LOGICAL hmsDriftCirc
      LOGICAL hmsDriftRect 
  
      REAL*8    uS(6)
      INTEGER i

      lost  = .TRUE.
      zz    = zT
 
      ! copy the target coordinates to the u vector 
      DO i=1,6
        u(i)=uT(i)
      ENDDO
      
      ! ----------------------------------------------------------- sive slit
      ! drift to sieve slit: z=1.27636m, 
      !   octagon edge m=2.545977, b=0.13502640m
      IF (hmsDriftOcta(u,zz,1.27636d00,0.0900176d00,0.0353568d00,
     >                                 2.5459770d00,0.1350264d00)) RETURN
  
      ! drift to back of sieve slit: z=1.33986m, dz=0.0635m
      !   octagon edge m=2.546569, b=0.141655189 m
      IF (hmsDriftOcta(u,zz,0.0635d00,0.0944368d00,0.0370839d00,
     >                             2.546569d00,0.141655189d00)) RETURN

      ! ------------------------------------------------------------------ Q1
      ! drift to mechanical entrance of Q1: z=1.4960m, dz=0.15614m
      IF (hmsDriftCirc(u,zz,0.15614d00,0.202575d00)) RETURN
      
      ! drift to Q1 entrance EFB: z=1.775805635m, dz=0.279805635m
      u(1)=u(1)+u(2)*0.279805635d00
      u(3)=u(3)+u(4)*0.279805635d00

      ! and save values
      DO i=1,6
         uS(i)=u(i)
      ENDDO       
      
      ! transport through Q1 fringe fields:
      CALL hmsApplyCOSY (u,1)
      IF (hmsCheckQuad(u,zz,0.279805635d00,0.202575d00)) RETURN
      
      ! transport through Q1, 1/5 at a time:
      DO i=1,4
        CALL hmsApplyCOSY (u,2)
        IF (hmsCheckQuad(u,zz,1.87838873d00/5.d00,0.202575d00)) RETURN
      ENDDO
      
      ! restore values:
      DO i=1,6
         u(i)=uS(i)
      ENDDO

      ! transport to Q1 exit EFB:
      CALL hmsApplyCOSY (u,3)
      IF (hmsCheckQuad(u,zz,1.87838873d00/5.d00,0.202575d00)) RETURN
      
      ! drift to Q1 mechanical exit: z=3.9340m, dz=0.279805635m
      IF (hmsDriftCirc(u,zz,0.279805635d00,0.202575d00)) RETURN

      ! ------------------------------------------------------------------ Q2
      ! drift to mechanical entrance of Q2: z=4.5610m, dz=0.6270m
      IF (hmsDriftCirc(u,zz,0.6270d00,0.29840d00)) RETURN
      
      ! drift to Q2 entrance EFB: z=4.887021890m, dz=0.326021890m
      u(1)=u(1)+u(2)*0.326021890d00
      u(3)=u(3)+u(4)*0.326021890d00
 
      ! and save values
      DO i=1,6
         uS(i)=u(i)
      ENDDO
      
      ! transport through Q2 fringe fields:
      CALL hmsApplyCOSY (u,4)
      IF (hmsCheckQuad(u,zz,0.326021890d00,0.29840d00)) RETURN
      
      ! transport through Q2, 1/5 at a time:
      DO i=1,4
        CALL hmsApplyCOSY (u,5)
        IF (hmsCheckQuad(u,zz,2.15595622d00/5.d00,0.29840d00)) RETURN
      ENDDO
      
      ! restore values:
      DO i=1,6
        u(i)=uS(i)
      ENDDO

      ! transport to Q2 exit EFB:
      CALL hmsApplyCOSY (u,6)
      IF (hmsCheckQuad(u,zz,2.15595622d00/5.d00,0.29840d00)) RETURN

      ! drift to Q2 mechanical exit: z=7.3690m, dz=0.326021890m
      IF (hmsDriftCirc(u,zz,0.326021890d00,0.29840d00)) RETURN
 
      ! ------------------------------------------------------------------ Q3
      ! drift to mechanical entrance of Q3: z=7.6610m, dz=0.2920m
      IF (hmsDriftCirc(u,zz,0.2920d00,0.29840d00)) RETURN
      
      ! drift to Q3 entrance EFB: z=7.990200290m, dz=0.329200290m
      u(1)=u(1)+u(2)*0.329200290d00
      u(3)=u(3)+u(4)*0.329200290d00
 
      ! save values
      DO i=1,6
         uS(i)=u(i)
      ENDDO
      
      ! transport through Q3 fringe fields:
      CALL hmsApplyCOSY (u,7)
      IF (hmsCheckQuad(u,zz,0.329200290d00,0.29840d00)) RETURN
      
      DO i=1,4
        CALL hmsApplyCOSY (u,8)
        IF (hmsCheckQuad(u,zz,2.14959942d00/5.d00,0.29840d00)) RETURN
      ENDDO
      
      ! and restore values:
      DO i=1,6
         u(i)=uS(i)
      ENDDO
      
      ! transport to Q3 exit EFB:
      CALL hmsApplyCOSY (u,9)
      IF (hmsCheckQuad(u,zz,2.14959942d00/5.d00,0.29840d00)) RETURN
     
      ! drift to Q3 mechanical exit: z=10.4690m, dz=0.329200290m
      IF (hmsDriftCirc(u,zz,0.329200290d00,0.29840d00)) RETURN
       
      ! -------------------------------------------------------------- Dipole
      ! drift to transition plate: z=11.058002m, dz=0.589002m
      IF (hmsDriftTPlate(u,zz,0.589002d00,0.30480d00,0.205232d00)) RETURN
      
      ! drift to opposite side of transition plate: z=11.092800m, dz=0.034798m
      IF (hmsDriftTPlate(u,zz,0.034798d00,0.30480d00,0.205232d00)) RETURN
      IF (hmsCheckDipole(u,zz,0.d00)) RETURN
      
      ! drift to D magnetic entrance: z=11.55m, dz=0.4572m
      u(1)=u(1)+u(2)*0.4572d00
      u(3)=u(3)+u(4)*0.4572d00
      IF (hmsCheckDipole(u,zz,0.4572d00)) RETURN

      ! save values:
      DO i=1,6
         uS(i)=u(i)
      ENDDO
      
      ! transport through 1/5 D with rotated entrance face:
      CALL hmsApplyCOSY (u,10)
      IF (hmsCheckDipole(u,zz,5.26053145d00/5.d00)) RETURN
      
      ! transport through 3/5 D with sector segments:
      DO i=1,3
        CALL hmsApplyCOSY (u,11)
        IF (hmsCheckDipole(u,zz,5.26053145d00/5.d00)) RETURN
      ENDDO
      
      ! restore values:
      DO i=1,6
         u(i)=uS(i)
      ENDDO

      ! transport through D (entrance to exit, fringe fields included):
      CALL hmsApplyCOSY (u,13)
      IF (hmsCheckDipole(u,zz,5.26053145d00/5.d00)) RETURN

      ! drift to transition plate: z=0.4572m, dz=0.4572m
      IF (hmsDriftTPlate(u,zz,0.457200d00,0.34290d00,0.205232d00)) RETURN
      IF (hmsCheckDipole(u,zz,0.d00)) RETURN
      
      ! drift to opposite side of transition plate: z=0.491998m,dz=0.034798m
      IF (hmsDriftTPlate(u,zz,0.034798d00,0.34290d00,0.205232d00)) RETURN
      
      ! drift to end of first piece of telescope: z=1.119378m, dz=0.62738m
      IF (hmsDriftCirc(u,zz,0.62738d00,0.338450d00)) RETURN
      
      ! drift to end of second piece of telescope: z=4.086098m, dz=2.96672m
      IF (hmsDriftCirc(u,zz,2.96672d00,0.384175d00)) RETURN
      
      ! drift to end of third piece of telescope: z=5.578398m, dz=1.4923m
      IF (hmsDriftCirc(u,zz,1.49230d00,0.460375d00)) RETURN

      ! ----------------------------------------------------------------- hut
      ! drift to focal plane.  This is the reference point for detector positions.
      u(1)=u(1)+u(2)*0.671602d00
      u(3)=u(3)+u(4)*0.671602d00
      zz = zz + 0.671602d00
      DO i=1,6
        uS(i)=u(i)
      ENDDO      
      
      ! drift to DC1 entrance: z=-0.51923-0.036=-0.55523m, dz=-0.55523m
      IF(hmsDriftRect(u,zz,-0.55523d00,-0.01670d00,0.565d00,-0.00343d00,0.26d00))RETURN
      ! drift to DC1 exit: z=-0.51923+0.054=-0.46523m, dz=0.090m
      IF(hmsDriftRect(u,zz, 0.09000d00,-0.01670d00,0.565d00,-0.00343d00,0.26d00))RETURN
      ! drift to DC2 entrance: z=0.29299-0.036=0.25699m, dz=0.72222m
      IF(hmsDriftRect(u,zz, 0.72222d00,-0.02758d00,0.565d00,-0.01653d00,0.26d00))RETURN
      ! drift to DC2 exit: z=0.29299+0.054=0.34699m, dz=0.090m
      IF(hmsDriftRect(u,zz, 0.09000d00,-0.02758d00,0.565d00,-0.01653d00,0.26d00))RETURN

      ! drift to S1X: z=0.7783m, dz=0.43131m
      IF (hmsDriftRect(u,zz,0.43131d00,0.015d00,0.6025d00,0.000d00,0.3775d00)) RETURN     
      ! drift to S1Y: z=0.9752m, dz=0.1969m
      IF (hmsDriftRect(u,zz,0.19690d00,0.000d00,0.6025d00,0.001d00,0.3775d00)) RETURN
      ! skip CK - no survey information
      ! drift to S2X: z=2.9882m, dz=2.013m
      IF (hmsDriftRect(u,zz,2.01300d00,0.004d00,0.6025d00,0.000d00,0.3775d00)) RETURN
      ! drift to S2Y: z=3.1851m, dz=0.1969m
      IF (hmsDriftRect(u,zz,0.19690d00,0.000d00,0.6025d00,0.013d00,0.3775d00)) RETURN

      ! drift to CAL: z=3.3869m, dz=0.2018m
      IF (hmsDriftRect(u,zz,0.20180d00,-0.134d00,0.6000d00,0.000d00,0.3000d00)) RETURN
       
      ! --------------------------------------------------------------- done 
      lost = .FALSE.
      DO i=1,6
        u(i)=uS(i)
      ENDDO      
       
      RETURN
      END

*------------------------------------------------------------------------
*------------------------------------------------------------------------
*------------------------------------------------------------------------
*
*       PART 2:  HMS Reconstruction (Backward Tracking; Focal Plane to Target)
*      -=======-       
*
*	Reconstruction (backward tracking) in the Jlab HMS hall C 
*       spectrometer using reconstruction and forward COSY matrices 
*       (including the effects of beam offsets (out-of plane)) 
*   
*       Both the normal in-plane scattering and the more special 
*       out-of-plane scattering are handeled. The later makes use 
*       of the forward COSY matrices. The algorithm was tested for
*       beam offsets in the range of cm (up or below the 
*       nominal scattering plane) 
* 
*       Supplies:
*         hmsInitRecon (map,p0) 
*           load the reconstruction maps
*         hmsInPlane (u,uT,ok)
*           reconstruction of the target coordinates 
*          (delta, dx/dz, y, dy/dz) at z=0
*         hmsOutOfPlane (u,x,uT,ok)
*           reconstruction of the target coordinates 
*           (delta, dx/dz, y, dy/dz) at z=0 including the 
*           vertical beam offset
*
*       Note: - Before calling hmsReconInPlane or hmsReconOutOfPlaneAccept 
*               the reconstruction map has to be loaded by a call to 
*               hmsInitRecon
*             - Before calling hmsReconOutOfPlane the forward transport 
*               maps have to be loaded by a call to hmsInitForward
*------------------------------------------------------------------------
 
      SUBROUTINE hmsInitRecon (filename) 
      IMPLICIT  NONE
      CHARACTER filename*(*)
      
* --  load the HMS reconstruction map  
*
*     Parameter:
*       filename    I : name of the reconstruction map
*       p0          I : momentum the spectrometer is set to
*
*     File format:
*      Ax' Ay Ay' Ad XxYy
*     
*       Ax'  : coeff's giving target x'
*       Ay   : coeff's giving target y
*       Ay'  : coeff's giving target y'
*       Ad   : coeff's giving the delta
*       X    : exponent of the focal plane x
*       x    : exponent of the focal plane x'
*       Y    : exponent of the focal plane y
*       y    : exponent of the focal plane y'

      ! matrix elements for the reconstruction
      INTEGER    NTERMS
      PARAMETER (NTERMS=1000)
 
      INTEGER   e1(NTERMS),e2(NTERMS),e3(NTERMS),e4(NTERMS),num,order
      REAL*8      c1(NTERMS),c2(NTERMS),c3(NTERMS),c4(NTERMS) 
      COMMON    /hmsRecon/num,order,e1,e2,e3,e4,c1,c2,c3,c4   
 
      ! other variables
      CHARACTER line*256
      INTEGER   i,eof,e5 
 
      ! read the necessary coeffs from the map data file
 10   FORMAT (A)
 20   FORMAT (1X,4G16.9,1X,5I1)
      
      OPEN (UNIT=97,FILE=filename,STATUS='OLD')

      num   = 1
      order = 0
              
      READ (97,10,IOSTAT=eof) line

      DO WHILE (eof .EQ. 0) 
        IF ((line(1:1) .NE. '!')     .AND. 
     >      (line(1:4) .NE. ' ---')  .AND.
     >      (line(1:2) .NE. 'h_'))   THEN
          ! read the reconstruction coefficents
          READ (line,20,IOSTAT=eof)  
     >      c2(num), c3(num), c4(num), c1(num),  
     >      e1(num), e2(num), e3(num), e4(num), e5  
          IF (((c1(num) .NE. 0.0d00) .OR. (c2(num) .NE. 0.0d00) .OR.
     >         (c3(num) .NE. 0.0d00) .OR. (c4(num) .NE. 0.0d00)) .AND.
     >        (e5 .EQ. 0) .AND. (eof .EQ. 0)) THEN
            i = e1(num) + e2(num) + e3(num) + e4(num) 
            IF (i .GT. order) order = i
            num = num+1
          ENDIF
        ENDIF
         
        IF (eof .EQ. 0) READ (97,10,IOSTAT=eof) line  
      ENDDO
      num = num-1
       
      CLOSE (97)
      
      RETURN
      END

*------------------------------------------------------------------------
      
      SUBROUTINE hmsReconOffset (uT,du)
      
      IMPLICIT NONE
      
      REAL*8 uT(6), du(4)
      
* --  calculates the focal plane correction for a given 
*     set of target coordinates
*
*     Parameter
*       uT   I : target coordinates (lab)
*                  uT(1,2) : x [m], dx/dz = out of plane coord. (downwards) 
*                  uT(3,4) : y [m], dy/dz = inplane coord. (perp. on x,z)
*                  uT(5)   : z [m] = in axis coordinate (towards HMS)
*                  uT(6)   : delta = relative deviation of the particle 
*                                    momentum from p0
*       du   O : correction values for the focal plane quantities
*                  du(1,2)  : correction in x [m], dx/dz 
*                  du(3,4)  : correction in y [m], dy/dz 
*
* frw 9/2000
* adjusted and repaired code in course of migration to g77
*
* some notes on how this routine works, since it does not obviously
* match the description given on MM's web page:
*
* Based on the web page documentation, this subroutine is supposed to transform
* the current guess at the reconstructed target coordinates back to the focal
* plane twice, once with x_target = beam value and once for x_target=0 and
* take the difference between the results as the correction du
* This code takes advantage of numerous cancellation and thereby reduces the 
* calculation effort:
* the sum below (DO i = first(20), last(20) ...) would be evaluated for both
* values, uT(1)=x and uT(1)=0.  The difference would be du.
* if e1(i)=0 then the result is independent of the value of uT(1) and since 
* that's the only difference between the 2 evaluations, these terms cancel
* in the difference giving du, thus they are omitted.
* similarly, if e1(i)<>0, the evaluation of the terms in the sum for uT(1)=0
* would result in all terms being 0, thus subtracting them from the sum for 
* uT(1)=x changes nothing.  Therefore, the difference between the transform 
* for uT(1)=x and uT(1)=0 is simply equal to the transform for uT(1)=x if the
* elements for which e1(i)=0 are skipped.
 
      include 'trans_map.inc'
  
      ! other variables
      REAL*8    a,uu1(0:10),uu2(0:10),uu3(0:10),uu4(0:10),uu6(0:10) 
      INTEGER i 
      
      uu1(0) = 1.d00
      uu2(0) = 1.d00
      uu3(0) = 1.d00
      uu4(0) = 1.d00
      uu6(0) = 1.d00
   
      ! drift backwards to z=0 
      uu1(1)=uT(1)-uT(2)*uT(5)
      uu3(1)=uT(3)-uT(4)*uT(5)
      
      uu2(1) = uT(2) 
      uu4(1) = uT(4)
      uu6(1) = uT(6)
     
     
      ! calculate the powers of the COSY coordinates
      DO i=1,4
        du(i) = 0.d00
      ENDDO
        
      DO i=2,order
        uu1(i)=uu1(i-1)*uu1(1)
        uu2(i)=uu2(i-1)*uu2(1)
        uu3(i)=uu3(i-1)*uu3(1)
        uu4(i)=uu4(i-1)*uu4(1)
        uu6(i)=uu6(i-1)*uu6(1)
      ENDDO
 
      ! calculate the focal plane offsets 
      
      DO i = first(20), last(20)

        IF (e1(i) .NE. 0) THEN
        
          a = uu1(e1(i)) * uu2(e2(i)) * uu3(e3(i)) 
     >                   * uu4(e4(i)) * uu6(e6(i))
     
          du(1) = du(1) + c1(i) * a
          du(2) = du(2) + c2(i) * a
          du(3) = du(3) + c3(i) * a
          du(4) = du(4) + c4(i) * a
          
        ENDIF

      ENDDO
                   
      RETURN
      END 

*------------------------------------------------------------------------

      SUBROUTINE hmsReconInPlane (u,uT,ok)
      IMPLICIT NONE
      REAL*8     u(4),uT(6)
      LOGICAL  ok
      
* --  performs the reconstruction of the target coordinates 
*     (delta, dx/dz, y, dy/dz) at z=0
*     
*     Parameter: 
*       u      I : focal plane coordinates (lab)  
*                    u(1,2)  : x [m], dx/dz = out of plane coord. (downwards) 
*                    u(3,4)  : y [m], dy/dz = inplane coord. (perp. on x,z)
*       uT     O : target coordinates (lab)
*                    uT(1,2) : x [m], dx/dz = out of plane coord. (downwards) 
*                    uT(3,4) : y [m], dy/dz = inplane coord. (perp. on x,z)
*                    uT(5)   : z [m] = in axis coordinate (towards HMS)  
*                    uT(6)   : delta (relative deviation of the particle 
*                                     momentum from p0)
*       ok   IO  : status variable 
*                   - if false no action is taken 
*                   - set to false when no reconstruction is found 

      ! matrix elemnts needed for calculating the focal plane offset 
      INTEGER    NTERMS
      PARAMETER (NTERMS=1000)
 
      INTEGER   e1(NTERMS),e2(NTERMS),e3(NTERMS),e4(NTERMS),num,order
      REAL*8      c1(NTERMS),c2(NTERMS),c3(NTERMS),c4(NTERMS)  
      COMMON    /hmsRecon/num,order,e1,e2,e3,e4,c1,c2,c3,c4   

      ! other variables
      REAL*8    a,uu1(0:10),uu2(0:10),uu3(0:10),uu4(0:10) 
      INTEGER i 


      DO i=1,6
        uT(i)  = 0.d00
      ENDDO
         
      ! calculate the powers of the focal plane coordinates
      
      uu1(0) = 1.d00
      uu2(0) = 1.d00
      uu3(0) = 1.d00
      uu4(0) = 1.d00
      
      uu1(1) = u(1)
      uu2(1) = u(2)
      uu3(1) = u(3)
      uu4(1) = u(4)
      DO i=2,order
        uu1(i)=uu1(i-1)*uu1(1)
        uu2(i)=uu2(i-1)*uu2(1)
        uu3(i)=uu3(i-1)*uu3(1)
        uu4(i)=uu4(i-1)*uu4(1)
      ENDDO
      ! calculate the target coordinates 
      DO i = 1, num
        a = uu1(e1(i)) * uu2(e2(i)) * uu3(e3(i)) * uu4(e4(i))
        
        uT(6) = uT(6) + c1(i) * a
        uT(2) = uT(2) + c2(i) * a
        uT(3) = uT(3) + c3(i) * a
        uT(4) = uT(4) + c4(i) * a 
      ENDDO

      ok = ((ABS(uT(2)) .LT. 1.) .AND. (ABS(uT(3)) .LT. 1.) .AND.
     >      (ABS(uT(4)) .LT. 1.) .AND. (ABS(uT(6)) .LT. 1.)) 
           
      RETURN
      END   
      
*------------------------------------------------------------------------

      SUBROUTINE hmsReconOutOfPlane (u,x,uT,ok)
      IMPLICIT NONE
      REAL*8     u(4),x,uT(6)
      LOGICAL  ok
      
* --  performs the reconstruction of the target coordinates 
*     (delta, dx/dz, y, dy/dz) for the out of plane case
*     
*     Parameter: 
*       u      I : focal plane coordinates (lab)  
*                    u(1,2)  : x [m], dx/dz = out of plane coord. (downwards) 
*                    u(3,4)  : y [m], dy/dz = inplane coord. (perp. on x,z)
*       x      I : x-offset at the target [m] (z=0) (out of plane; downwards)
*       uT     O : target coordinates (lab)
*                    uT(1,2) : x [m], dx/dz = out of plane coord. (downwards) 
*                    uT(3,4) : y [m], dy/dz = inplane coord. (perp. on x,z)
*                    uT(5)   : z [m] = in axis coordinate (towards HMS)  
*                    uT(6)   : delta (relative deviation of the particle 
*       ok   IO  : status variable 
*                   - if false no action is taken 
*                   - set to false when no reconstruction is found 
*                                     momentum from p0)
      
      REAL*8       eps           
      PARAMETER (eps = 0.0005d00) ! accuracy in delta
      
      REAL*8       dd,du(4),u0(4) 
      INTEGER    n
             

      CALL hmsReconInPlane (u,uT,ok)     ! first guess
      uT(1) = x
c
c
      du(1) = 0.d00
      du(2) = 0.d00
      du(3) = 0.d00
      du(4) = 0.d00
           
      dd = 1.d00
      n  = 0
      DO WHILE ((ABS(dd) .GT. eps) .AND. (n .LT. 10) .AND. ok)
      
        CALL  hmsReconOffset (uT,du)

        u0(1) = u(1)-du(1)
        u0(2) = u(2)-du(2)  
        u0(3) = u(3)-du(3)  
        u0(4) = u(4)-du(4)  
 
        dd = uT(6)
        
	CALL hmsReconInPlane (u0,uT,ok)
        uT(1) = x

        dd = dd-uT(6)
        n  = n+1
        
      ENDDO 

      IF (ABS(dd) .GT. eps) ok = .FALSE.    ! not converged

      RETURN
      END


