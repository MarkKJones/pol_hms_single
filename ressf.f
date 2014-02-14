CCC  Subroutine to get F1 and FL from fits to L/T cross sections.    CCC
CCC  The subroutine resmod.f is required.                            CCC

cc   Modified from the original routine (done by M.Jones, S.tajima)
cc
cc   Variable names changed to avoid confusion (F_one, F_L)
cc   Initialize variables (F_one, R, F2)
cc   Factor for F_one and F_L changed to  1/0.389e3
cc   Add lines for calculating F2 and R, and return these instead of F1 and FL
cc   Set R=F2=0 if F_one=0
cc

      SUBROUTINE ressf(Q2,W2,F2,R)

      IMPLICIT NONE

      real*8 w2,q2,xval1(50),xvall(50),temp(4)
      real*8 mp,mp2,pi,alpha,xb,F_one,F_L,F2,r,g2
      integer i,lu
      logical first/.true./
      integer ntp_events
      common /ntpevents/ ntp_events
 
      F_one=0.d0
      F_L=0.d0
      R=0.d0
      F2=0.d0
c
      mp = .9382727
      mp2 = mp*mp
      pi = 3.141593
      alpha = 1./137.036

      lu=63
      if(first) then
        first=.false.
        open(unit=lu,file='f1parms.dat',status='old')
        do i=1,50
          read(lu,*) temp(1)  ! par #
          read(lu,*) xval1(i) ! starting value
          read(lu,*) temp(2)   ! initial step (0 means fixed parm)
          read(lu,*) temp(3) ! low limit
          read(lu,*) temp(4) ! high limit 
        enddo
        close(lu)

        open(unit=lu,file='flparms.dat',status='old')
        do i=1,50
          read(lu,*) temp(1)  ! par #
          read(lu,*) xvall(i) ! starting value
          read(lu,*) temp(2)   ! initial step (0 means fixed parm)
          read(lu,*) temp(3) ! low limit
          read(lu,*) temp(4) ! high limit 
        enddo
        close(lu)
      endif

      xb = q2/(w2+q2-mp2)
      if ( w2 .ge. (0.9382727+0.136)**2) then
       call resmod(1,w2,q2,xval1,F_one)   !!! Transverse Cross Section    !!!

       call resmod(2,w2,q2,xvall,F_L)   !!! Longitudinal Cross Section  !!!
      endif
cc      if (F_one .ne. F_one .or. r.ne.r .or. F2.ne.F2 .or. F_L.ne.F_L) then
cc        print *, 'ressf.f (just after resmod): F_one=',F_one, ' F_L=',F_L
cc      endif

CCC  Now convert to structure functions  CCC

      if ( F_one .eq. 0.d0 ) then
        R=0.d0  
        F2=0.d0
      else
        F_one = F_one*(w2-0.9382727*0.9382727)/8./pi/pi/alpha/0.389e3
        F_L = F_L*2*xb*(w2-0.9382727*0.9382727)/8./pi/pi/alpha/0.389e3
        R = F_L/ (2.D0 * XB * F_one)
	g2 = (2.d0*xb*mp)**2/q2
	F2 = (F_L+2.d0*xb*F_one)/(1+g2) 
      endif

cc      if (F_one .ne. F_one .or. r.ne.r .or. F2.ne.F2 .or. F_L.ne.F_L) then
cc      print *,'ressf.f: F_one=',F_one,' F2=',F2,' F_L=',F_L,' r=',r,' w2=',w2
cc             ALL of these became NaN !!!!
cc      endif
    
      end
