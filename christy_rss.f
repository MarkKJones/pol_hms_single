!      subroutine christy(W2,Q2,F1,R)
      subroutine christy_rss(q2,w2,f2,r)
!------------------------------------------------------------------------
! subroutine to return proton structure function F1 and ratio R=sigl/sigt
!
! inputs are electron missing mass squared W2 (GeV**2) 
!            momentum trasfer squared Q2 (GeV**2)
! inputs and outputs are Real*8
! the file christy.dat is needed to use this subroutine
! Note: if W2<1.155 GeV**2, values of zero are returned (below threshold)
! Fit done by Eric Christy 11/04
! Reference this fit as ???       
!------------------------------------------------------------------------
      IMPLICIT NONE

      real*8 w2,q2,xval1(41),xvall(41),temp(4)
      real*8 mp,mp2,pi,alpha,xb,F1,FL,R
      real*8 f2,g2
      integer i
      integer lu
      logical first/.true./
      
      mp = .9382727d0
      mp2 = mp*mp
      pi = 3.141593d0
!      alpha = 1./137.
      alpha = 1.d0/137.04d0
	lu=63
      if(first) then
        first=.false.
!        open(unit=15,file='christy.dat',status='old')
        open(unit=lu,file='christy.dat',status='old')
        do i=1,41
          read(lu,*) temp(1)  ! par #
          read(lu,*) XVAL1(i) ! starting value
          read(lu,*) temp(2)   ! initial step (0 means fixed parm)
          read(lu,*) temp(3) ! low limit
          read(lu,*) temp(4) ! high limit 
        enddo
        do i=1,41
          read(lu,*) temp(1)  ! par #
          read(lu,*) XVALL(i) ! starting value
          read(lu,*) temp(2)   ! initial step (0 means fixed parm)
          read(lu,*) temp(3) ! low limit
          read(lu,*) temp(4) ! high limit 
        enddo
        close(lu)
      endif

      F1=0.d0
      R=0.d0
      f2=0.d0
      if(w2.lt.1.155) return

      xb = q2/(w2+q2-mp2)
      if(xb.le.0.d0) return

      call resmod2(1,w2,q2,XVAL1,F1)
      call resmod2(2,w2,q2,XVALL,FL)

       if(F1.le.0.0) return

      F1 = F1/8.d0/pi/pi/alpha/1000.d0
      F1 = F1*(w2-mp2)
      FL = FL/8.d0/pi/pi/alpha/1000.d0
      FL = FL*(w2-mp2)*2.d0*xb
      R = FL/ (2.D0 * XB * F1)
	g2 = (2.d0*xb*mp)**2/q2
	f2 = (fl+2.d0*xb*f1)/(1+g2) 
      return
      end

      SUBROUTINE RESMOD2(sf,w2,q2,xval,fn) 

      IMPLICIT NONE
      REAL*8 W,w2,q2,mp,mp2,xb,fn,xval(41),mass(4),width(4)
      REAL*8 height(4),fn_del,fn_s11,fn_f15,rescoef(4,3)
      REAL*8 nr_coef(4,4),wdif,fn_nr,fn_4,w2temp,wtemp
      REAL*8 roper_mass,roper_width,roper_height
      REAL*8 roper_mparm,roper_exp,mq2(4)
      INTEGER i,j,num,sf


      mp = 0.9382727
      mp2 = mp*mp
      W = sqrt(w2)
      w2temp = w2 - .9382727*.9382727
      wtemp = sqrt(w2temp)
      wdif = w - (0.937272 + 0.137)
      xb = q2/(q2+w2-mp2)
      if(sf.EQ.1) xb = 1.d0
      fn_nr = 0.0d0
      do i=1,4
        height(i) = 0.
      enddo
      
      do i=1,4
        mass(i) = xval(i)
      enddo

      do i=1,4
        mq2(i) = xval(36 + i)
      enddo

      if(q2.GT.0.1) then

        mass(2) = mass(2)*exp(-(q2-0.1)/mq2(2)) 
     &            + mq2(1)*(1.-exp(-(q2-0.1)/mq2(2))) 

        mass(3) = mass(3)*exp(-(q2-0.1)/mq2(4)) 
     &            + mq2(3)*(1.-exp(-(q2-0.1)/mq2(4)))

      endif   
      do i=1,4
        width(i) = xval(4+i)
      enddo
      num = 0
      do i=1,4
       do j=1,3
         num = num + 1
         rescoef(i,j)=xval(8 + num)
c           write(6,*) i,j,num,rescoef(i,j)
c           height(i) = height(i)+rescoef(i,j)*q2**(-1.*(float(j-1)))
         enddo
         height(i) = rescoef(i,1)*
     &             (1.+q2/rescoef(i,2))**(-1.*rescoef(i,3))
c         if(w2.LT.1.35) height(i) = height(i)/(w2-1.30)    
      enddo
      num = 0     
      do i=1,4
       do j=1,4
         num = num + 1
         nr_coef(i,j)=xval(20 + num)
c         write(6,*) i,j,num,nr_coef(i,j)         
       enddo
      enddo

      do i=1,5
       roper_mass = xval(37)
       roper_width = xval(38)
       roper_height = xval(39)
       roper_mparm = xval(40)
       roper_exp = xval(41)
      enddo
c      write(6,*) "constant coef's are:  ",height

CC   Calculate Breit-Wigners for the 3 resonance regions   CC

      fn_del = width(1)/((W-mass(1))**2 
     &               + 0.25*width(1)*width(1))
      fn_s11 = width(2)/((W-mass(2))**2 
     &               + 0.25*width(2)*width(2))
      fn_f15 = width(3)/((W-mass(3))**2 
     &               + 0.25*width(3)*width(3))
      fn_4   = width(4)/((W-mass(4))**2 
     &               + 0.25*width(4)*width(4))

      fn_del = height(1)*fn_del
      fn_s11 = height(2)*fn_s11
      fn_f15 = height(3)*fn_f15
      fn_4   = height(4)*fn_4

c      roper_height = roper_height*(1.+q2/roper_mparm)**(-1.*roper_exp)
c      fn_roper = roper_width/((W-roper_mass)**2 
c     &               + 0.25*roper_width*roper_width) 
c      fn_roper = roper_height*fn_roper


      do i=1,4
       do j=1,4
c         fn_nr = fn_nr+
c     &      nr_coef(i,j)*q2**(float(j-1))*sqrt(wdif**(float(i)))

         fn_nr = fn_nr+
     &      nr_coef(i,j)*q2**(float(j-1))*sqrt(wdif**(float(i)))
     &                  /w2temp/xb

c         write(6,*) sf

c         if(sf.EQ.2) fn_nr = fn_nr/xb
                         
       enddo
      enddo

      fn = fn_del + fn_s11 + fn_f15 + fn_4 + fn_nr
      if(sf.EQ.2) then
        fn = fn*(1.-exp(-q2/roper_mass))
      endif

      !write(6,*) "IN model:  ",w2,q2,fn
!      write(69,*) w2,q2
      RETURN 
      END 


