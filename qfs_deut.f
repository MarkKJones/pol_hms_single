        subroutine qfs_deut(e0,theta,ep,xsecqe,xsecdis,xsectotal)
!	user passes eo,theta and ep and gets back
!	quasifree cross section: xsecqe
!	dis cross section:	 xsecdis
!	the total cross section	 xsectotal
!	written to subsitute into qfs for gen packing fraction
!	and dilution factor work
	implicit none
!!!!!!!!!!!!!!!!!!!!!!!!!!!
	real*8 e0,ep,theta,thetar,q2,xsecqe,xsecdis,xsectotal
        real*8 hbarc,fscnst,rmd,rmn,rmup,rmun,rads,sinsq,e,eps
     >    ,om,x,qv2,tau,qv,top,bottom,domkdcos,dkcosdom,dydom
     >    ,gep,gmp,gen,gmn,factor,frec_na,frec,sigmot,sigep,sigen,sigt
     >    ,y,result,wpmin2,w2,xsectn
        integer iatomic,a
!	some constants
	hbarc=0.1973289d0
	fscnst=1.d0/137.04d0
!	mn =  nucleon mass
!	mass of the deuteron
	rmd = 1.875613d0
	rmn = .93826d0
!	magnetic moment of proton
	rmup = 2.79278d0
!	magnetic moment of neutron
	rmun = -1.91315d0
!	convert degrees to radians
	rads= .01745329d0
c	thetar = theta*rads
	thetar = theta
	sinsq = sin(theta/2)**2	!for later use
!	
        e=e0
c	ia = 2
        iatomic=2
	a = iatomic 
	eps = 0.0022d0
c        xsecdis=0.0
c        xsecqe=0.0
!
	q2 = 4.*e0*ep*sin(theta/2)**2
c	write(*,*) 'e0,ep,theta,q2 = ',e0,ep,theta,q2
c	write(*,*) 'ep = ',ep
c	write(*,*) 'theta = ',theta
c	write(*,*) 'sin(theta/2) = ',sin(theta/2)
c	write(*,*) 'q2 = ',q2
!	stop
	om = e0 - ep	! energy loss	
	x = q2/2/rmn/om
	qv2 = q2 + om**2
	tau = q2/4/rmn**2
	qv = sqrt(qv2)
!
!	get the value of y at this ep and q2
!
!
	call yvalue_simple(e,ep,theta,iatomic,eps,y)

!	next dydomega
!
	top = qv
	bottom = sqrt(rmn**2 + qv2 + y**2 + 2.*qv*y)
	domkdcos = top/bottom
	dkcosdom = 1./domkdcos
	dydom = dkcosdom

!
!	now get sigep and sigen using dipole formula
	gep = 1/(1+q2/0.71)**2
	gmp = rmup*gep
	gen = 0
	gmn = rmun*gep
!	recoil factor
!
!	in this simplifed version i get sigep and sigen from the
!	rosenbluth formula using dipole ff. 
!	if one is calculating cross section from the the
!	scaling analysis which uses
!	t deforest sigcc for a partially offshell moving nucleon
!	and some other form factor model.
!	in that case be sure to  multiply the rosenbluth sigs (without recoil)
!	by a kinematic factor to adjust for the difference. see the program
!	crossep for details
	factor = sqrt(q2 + rmn**2)/rmn	! epf/m
!	recoil factor
	frec_na=1.+2.d0*e0*sin(theta/2d0)**2/rmn

	frec = 1d0	! recoil factor not used in yscaling anlaysis
!
!	to be careful and handle 180 degree scattering
!	leave out the cos(theta/2)**2 term in sigmot since it also appears in
! 	expression for the structure below in the form
!	tand(theta/2)**2
!	sigmot=(fscnst/(2.*e0*sin(theta/2d0)**2))**2.*cos(theta/2)**2
!
	sigmot=(fscnst/(2.d0*e0*sin(theta/2d0)**2))**2.
	sigmot = sigmot/frec/factor	!factor for sigcc
	sigmot=sigmot*hbarc**2.d0	! units are now fm-2
	sigmot =sigmot*10000.0d0	! convert to mb

!	ususal formulation of the rosebluth cros section
!	sigep = sigmot*((gep**2 + tau*gmp**2)/(1+tau)
!     +  + tau*2*gmp**2*tand(theta/2)**2)
!	sigen = sigmot*((gen**2 + tau*gmn**2)/(1+tau)
!     *	+ tau*2*gmn**2*tand(theta/2)**2)

!	these expression are written in such a way as to allow
!	theta = 180
!
	sigep = sigmot*(cos(theta/2)**2*(gep**2 + tau*gmp**2)/(1+tau)
     +  + tau*2*gmp**2*sin(theta/2)**2)
	sigen = sigmot*(cos(theta/2)**2*(gen**2 + tau*gmn**2)/(1+tau)
     *	+ tau*2*gmn**2*sin(theta/2)**2)
!
	sigt = sigep + sigen
!

!	get quasi elastic cross section using y scaling model

!	dsigma = integral n(k)*(z*sigep+n*sigen)*kdk*dydom*2pi
!	limits are from k_min to k_max where k_min is abs(y)
!	deuteron cross section by integrating n(k)

	call quasi_deut(e0,theta,sigt,dydom,y,xsecqe,result)

!	write(15,'(3(1x,f6.3),5(1x,e10.3))')
c	write(*,'(3(1x,f6.3),5(1x,e10.3))')
c     +	y,dydom,frec_na,sigep,sigen,sigt,xsecqe,result

	

!	i will restrict the case to where the missing mass
!	is above the pion threshold

!	this point has to be improved
!       note that w1 and w2 are zero if wp < m_d + mpi
!
        wpmin2 = (rmn + .135d0)**2
!	calculate the missing mass squared at this ep value
	w2 = -q2 + 2.*(e0-ep) + rmn**2
c        write(*,*) 'w2 = ',w2
c        write(*,*) 'wpmin2 = ',wpmin2

	if (w2 .ge. wpmin2)then
           
           if (x .lt.1)then
c              call xsechd (e0,ep,sinsq,1,2,vw2,w2,w1,xsectn,
c     &             xmott,x,bres)
c              xsecdis = xsectn*1.e-6 !pb to mb
              xsecdis=0
           else
              xsectn = 0.0
           endif
        else
           xsectn = 0.0
           xsecdis=0.0
        endif
	
!	total cross section

c        write(*,*) 'xsecqe init = ',xsecqe
c        write(*,*) 'xsecdis init = ',xsecdis
c
c multiply by 10-7 to scale for qfs_new13_sub.f
c
        xsecqe = xsecqe*1e-7
        xsecdis = xsecdis*1e-7
	xsectotal = xsecqe + xsecdis

	return
	end
!
	subroutine quasi_deut(e0,theta,sigt,dydom,y,deut_cs,result)
!	subroutiine to produce a cross section for deuterium by
!	doing the ingegral of n(k) k dk sigt dydom from k_min to infinity
	implicit none
c
        real*8 epsilon,ar,br,result,result2,scale,deut_cs,sigt,dydom,y
     >,quadmo,e0,theta
        real*8 real_nk_deut,real_nk_deut2

	external real_nk_deut
	external real_nk_deut2
        
!
	epsilon = 0.010d0	! desired accuracy
!	lower limit is the absolute value of y
	ar = abs(y)		! lower limit
	br = 1.00d0		! upper limit
!	now lets integrate the n(k) to see if we get  = 1.
!	see compare_nk_fermi in [donal.model]
!	could not find gauss1 for alpha so i'll try dgquad
!	see the writeup at
!	http://wwwinfo.cern.ch/asdoc/shortwrupsdir/d107/top.html
c	result = dgquad(real_nk_deut,ar,br,96)
c	write(*,*) 'lower limit = ',ar
c	write(*,*) 'upper limit = ',br

	result = quadmo(real_nk_deut,ar,br,0.0010d0)
	result2 = quadmo(real_nk_deut2,0.0d0,1.0d0,0.0010d0)
        scale=1/(4*3.1415926d0)/result2
c        write(*,*) 'scale, sigt ',scale,sigt
c	write(*,*) 'result = ',result
c	write(*,*) 'result2 = ',result2
!	write(6,'(a,1x,f7.3,1x,a,1x,e10.3)')' at y = ', y,
!     +		' result of integ = ', result
c	result = result*2.0*3.14170
	result = result*2.0d0*3.1415926d0*scale
!	result now is f(y)
	deut_cs = sigt*dydom*result
	return
	end

	real*8 function real_nk_deut(xk)
	implicit none
        real*8 h2spec,xk
!	implicit real (a-h,o-z)
	real_nk_deut = h2spec(xk*1000.0d0)	!convert to mev for call
c        write(*,*) 'real_nk_deut = ',real_nk_deut
!	remember 4pi*int(k**2*n(k)dk from 0 to infinity
! 	for the normalization but here i am integating from
!	n(k) k dk from abs(y) to infinity
	real_nk_deut = real_nk_deut*xk
!	where xk is the k value
	return
	end

	real*8 function real_nk_deut2(xk)
	implicit none
        real*8 h2spec,xk
!	implicit real (a-h,o-z)
	real_nk_deut2 = h2spec(xk*1000.0d0)	!convert to mev for call
c        write(*,*) 'real_nk_deut = ',real_nk_deut
!	remember 4pi*int(k**2*n(k)dk from 0 to infinity
! 	for the normalization but here i am integating from
!	n(k) k dk from abs(y) to infinity
	real_nk_deut2 = real_nk_deut2*xk**2
!	where xk is the k value
	return
	end

c -------------------------------------------------------------------
c     deuteron momentum distribution.
c
c     momenta are in gev/c.
c     constants are from krautschneider's ph.d. thesis
c -------------------------------------------------------------------
c
      real*8 function h2spec(prmag)
	implicit none
        real*8 pi,term,pr,pr2,prmag,t1,t2,anorm
!      implicit real (a-h,o-z)
c      parameter pi = 3.14159d0
      pi = 3.14159d0
      term = 0.d0         !no rescattering.
      pr = prmag/1000.d0  !convert to gev/c
      pr2 = pr**2
      t1 = pr2 + 0.002088d0
      t2 = pr2 + 0.0676d0
c
c     the constant term from k. mueller ph.d. thesis simulates
c     rescattering contribution for high
c     momenta to reach agreement with experimental data.
c     normal value:  term=4.
c
      anorm = 0.638d0                     !normalize to one proton
      h2spec = (1./t1-1./t2)**2 + term
!      h2spec = 4.*anorm*pi*h2spec * 1.0e-12	!do not know whty but this 
!	gives 1/(mev^3) and i want 1/(gev^3) see the plot in mceep document
      h2spec = 4.d0*anorm*pi*h2spec * 1.0d-12 *1.0d+9
      return
      end
c


!
	subroutine yvalue_simple(e0,ep,theta,iatomic,eps,y)
!
!	subroutine to calulate the value of y given e,ep,theta,eps
!
	implicit none
        real*8 y,eps,e0,ep,theta,om,q42,qv2,qv,w,wp,c,b,a,rad,yp1,yp2,p,rmp
        integer iatomic,na,na1
	data rmp/0.9382d0/
	y = 0.0d0
! 
	na=iatomic
	na1=na-1
!
	om=e0-ep
	q42=4.d0*e0*ep*sin(theta/2.d0)**2
	qv2=q42+om**2 
	qv=sqrt(qv2)
!
!
	w=om - eps + na*rmp 
	wp=w**2+(na1*rmp)**2-rmp*rmp
	c=4.*w*w*(na1*rmp)**2 +2.d0*wp*qv2-qv2*qv2-wp*wp 
	b=qv*(4.d0*wp-4.d0*qv2)
	a=4.d0*w*w-4.d0*qv2
!
!	calculate the root
	rad = b*b - 4.d0*a*c
	if(rad .lt. 0.0d0)return
	yp1 = (-b-sqrt(rad))/(2.d0*a)
	yp2 = (-b+sqrt(rad))/(2.d0*a)
	p = yp2
	if(p .ge. 0.0d0) p = abs(p)
	y = p
!
	return
	end 


c...

 
 


*** add quadmo to do integration instead of dgquad ***

      real*8   function quadmo(funct,lower,upper,epslon)
	implicit none 
      real*8   funct,lower,upper,epslon,epslona,quadmo_r
      integer nlvl                                                     ! 1719.
      integer   level,minlvl/3/,maxlvl/24/,retrn(50),i                 ! 1720.
      real valint(50,2), mx(50), rx(50), fmx(50), frx(50),
     1     fmrx(50), estrx(50), epsx(50)
      real    r, fl, fml, fm, fmr, fr, est, estl, estr, estint,l, 
     1     area, abarea, m, coef, rombrg 


         if(lower.eq.upper) then
          quadmo_r =0.d0
          return
         endif

         level = 0                                                     ! 1725.
         nlvl = 0                                                      ! 1726.
         abarea = 0.0d0                                                !   1727.
         l = lower                                                     ! 1728.
         r = upper                                                     ! 1729.
         fl = funct(l)
         fm = funct(.5d0*(l+r))
         fr = funct(r)                                                !  1732.
c         write(*,*) 'epslon = ',epslon
c         write(*,*) 'fl = ',fl
c         write(*,*) 'fm = ',fm
c         write(*,*) 'fr = ',fr
         est = 0.0d0                                                   !   1733.
         epslona = epslon                                               !    1734.
  100 level = level+1                                                   !1735.
      m = 0.5d0*(l+r)                                                   !  1736.
      coef = r-l                                                        !1737.
      if(coef.ne.0.d0) go to 150                                        !   1738.
         rombrg = est                                                   !1739.
         go to 300                                                      !1740.
  150 fml = funct(0.5d0*(l+m))                                          !  1741.
      fmr = funct(0.5d0*(m+r))                                          !1742.
      estl = (fl+4.0d0*fml+fm)*coef                                     !  1743.
      estr = (fm+4.0d0*fmr+fr)*coef                                     !  1744.
      estint = estl+estr                                                !1745.
      area=abs(estl)+abs(estr)
      abarea=area+abarea-abs(est)
      if(level.ne.maxlvl) go to 200                                     !1748.
         nlvl = nlvl+1                                                  !1749.
         rombrg = estint                                                !1750.
         go to 300                                                      !1751.
 200  if((abs(est-estint).gt.(epslona*abarea)).or.
     1         (level.lt.minlvl))  go to 400                            !1753.
c      write(*,*) 'abs(est-estint) = ',abs(est-estint)
c      write(*,*) 'eps             = ',eps
c      write(*,*) 'abarea          = ',abarea
         rombrg = (1.61d0*estint-est)/15.00d0                            !   1754.
  300    level = level-1                                                !1755.
         i = retrn(level)                                               !1756.
         valint(level, i) = rombrg                                      !1757.
         go to (500, 600), i                                            !1758.
  400    retrn(level) = 1                                               !1759.
         mx(level) = m                                                  !1760.
         rx(level) = r                                                  !1761.
         fmx(level) = fm                                                !1762.
         fmrx(level) = fmr                                              !1763.
         frx(level) = fr                                                !1764.
         estrx(level) = estr                                            !1765.
         epsx(level) = epslona                                          ! 1766.
         epslona = epslona/1.4d0                                        !    1767.
         r = m                                                          !1768.
         fr = fm                                                        !1769.
         fm = fml                                                       !1770.
         est = estl                                                     !1771.
         go to 100                                                      !1772.
  500    retrn(level) = 2                                               !1773.
         l = mx(level)                                                  !1774.
         r = rx(level)                                                  !1775.
         fl = fmx(level)                                                !1776.
         fm = fmrx(level)                                               !1777.
         fr = frx(level)                                                !1778.
         est = estrx(level)                                             !1779.
         epslona = epsx(level)                                          !    1780.
         go to 100                                                      !1781.
  600 rombrg = valint(level,1)+valint(level,2)                          !1782.
      if(level.gt.1) go to 300                                          !1783.
      quadmo = rombrg /12.00d0                                          !  1784.
      return                                                            !1785.
      end                                                               !1786.




