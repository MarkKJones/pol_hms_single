*qfs
      subroutine qfs(z,a,e,th,w,delnu,flag,sigma_qfs,elflag,scalein)
c
      implicit real*8 (a-h,o-z)
      character*1 dorg
      logical*4 flag
      real*8 elflag
      common /par/ eps,epsd,pf,spence
      common /switch/ dorg
      parameter(scale=1.d-26)
      parameter (wpi = 1.078d0)

      pm=939.0d0
      dm=1219.d0
      alph=1.d0/137.03604d0
      pi=acos(-1.d0)
c
      IF(A.EQ.1.) THEN
        PF=0.    
        EPS=0.
        EPSD=13.0
        dorg='y'
      ELSEIF(A.eq.2.) THEN
!        PF=77. 
!        EPS=2.23
!        PF=74.	!	OR 2/04
        PF=100.	!	OR 2/05
        EPS=2.22
        EPSD=13.0
      ELSEIF(A.EQ.4.) THEN
!        PF=180.  !152   
!        EPS=20.  !8
        PF=162.  ! 	OR 3/05   
        EPS=20.6  !	OR 3/05
        EPSD=13.0
      ELSEIF(A.EQ.9.) THEN
        PF=200.       
        EPS=20.
        EPSD=13.0
      ELSEIF(A.EQ.12.) THEN
        PF=230.    !205, 215, 220, 230=5/03     
        EPS=25. 
        EPSD=13.
      ELSEIF(A.EQ.14 .OR. A.EQ.15) THEN
        PF=225.      
        EPS=25.
        EPSD=13.
      ELSEIF(A.EQ.27.) THEN
        PF=235.       
        EPS=32.       
        EPSD=13.
      ELSEIF(A.EQ.59.) THEN
        PF=260.               
        EPS=35.
        EPSD=13.
      ELSEIF(A.EQ.64.) THEN
        PF=260.                
        EPS=35.
        EPSD=13.
      ELSEIF(A.EQ.184.) THEN
        PF=265.                
        EPS=43.
        EPSD=13.
      ELSE
        WRITE(*,*) 'pf, eps, epsd not specified for this material, stop'
        STOP
      ENDIF
c
      sigqfza=sigqfs(e,th,w,z,a,eps,pf,delnu)*scale
      
!	6/03 Use Hall C F2 R for proton inelastic sigma

      thr = th*pi/180.0d0	! 6/03
      qsq = 4.0d0*e*(e-w)*sin(thr/2.0d0)**2	! 6/03
      wsq = pm*pm + 2.d0*pm*w - qsq	! 6/03

      if((a.eq.1.and.z.eq.1).and.(wsq/1d6.le.3.9d0).and.elflag .eq. 0.d0) then	! 6/03

	  call hcf2r(qsq/1.d6,wsq/1.d6,f2,r)	! 6/03
	  siginel=sigmot(e,thr)*f2*(1.d0/w+
	1	2*w*(1+qsq/w**2)*(tan(thr/2))**2/(qsq*(1+r)))	! 6/03
	   siginel=siginel*scale	! 6/03

      else				! 6/03 A >1

        if(a.eq.2.and.z.eq.1) then	! 2/04

!	  print *,'e ',e,enu,ww,qsq,th
		call hallc2h(e,th,w,sigdelc,sigr1c,sigr2c,sigxc)	! 2/04
		sigda = sigdelc*scale	!	"
		sigr1a = sigr1c*scale	!	"
		sigr2a = sigr2c*scale	!	"
		sigxa = sigxc*scale		!	"
	       siginel = sigda+sigxa+sigr1a+sigr2a	!	"

!	    if (w.lt.wpi) then
!	    if (w.lt.1140) then	!	
!		sig2na=sig2n(e,th,enu,z,a,pf)*scale*0.7d0
!		sig2na=sig2n(e,th,enu,z,a,pf)*scale*1.4d0	!	"
!             sig2na=sig2n(e,th,w,z,a,pf)*scale*0.7
!	    else	!	"
		sig2na = 0.d0	!	"
!	    end if ! No 2N sig if W > W_pi_thr. 2/04

	 else	!	A > 2

         sigda=sigdel(e,th,w,a,epsd,pf,z)*scale
         sigxa=sigx(e,th,w,a)*scale
         sigr1a=sigr1(e,th,w,a,pf,z)*scale
         sigr2a=sigr2(e,th,w,a,pf,z)*scale
!         sig2na=sig2n(e,th,w,z,a,pf)*scale*0.7
!         siginel=sigda+sigxa+sigr1a+sigr2a+sig2na	! 6/03
	  siginel=sigda+sigxa+sigr1a+sigr2a
          siginel = siginel*scalein
!	  sig2na=sig2n(e,th,enu,z,a,pf)*scale*0.7
         sig2na=sig2n(e,th,w,z,a,pf)*scale*1.0d0		! 11/04

	 end if ! hallc 2h

!	   sig2na=sig2n(e,th,enu,z,a,pf)*scale*0.7
	  siginel=siginel+sig2na
c  
 
	end if	! hcf2r	2/04
!      end if						! 6/03
      
!      sig=sigqfza+sigda+sigxa+sigr1a+sigr2a+sig2na
	sig = sigqfza + siginel				! 6/03
c
        if ( (a.eq.1.and.z.eq.1) ) then
         if (  elflag .eq. 1) then
           sig = sigqfza
         else
           sig = siginel 
         endif
        endif
c	  
      sigma_qfs=sig        ! nonradiated qfs cross section
c
c      thcr=th*pi/180.
c      qms=4.*e*(e-w)*sin(thr/2.)**2
c      qvs=qms+w**2
c      ekappa=w-qms/2./pm
c      if(ekappa.gt.-pm/2.)then
c        cmtot=sqrt(pm**2+2.*pm*ekappa)
c      else
c        cmtot=pm
c      endif
c      flux=(alph/2./pi**2)*((e-w)/e)*((2.*pm*w-qms)/2./pm/qms)
c      polari=1./(1.+2.*qvs*tan(thr/2.)**2/qms)
c      flux=flux/(1.-polari)
c      if(ekappa.lt.0.) photsig=0.
c      photsig=sig/flux
c radiate
      sigrad=0.
      if(sig.gt.0.0d0.and.flag) then
	call radiate(e,th,w,z,a,sig,sigrad)
        sigma_qfs=sigrad
      end if
c
      return
      end

*fd
      real*8 function fd(qms,a)
      implicit real*8 (a-h,o-z)
      fd=1./(1.+qms/a**2)**2
      return
      end

*fm
      real*8 function fm(qms,a)
      implicit real*8 (a-h,o-z)
      fm=1./(1.+qms/a**2)
      return
      end

*fphenom
      real*8 function fphenom(qms)
      implicit real*8 (a-h,o-z)
      a1=.55d0
      a2=20.d0/1.d6
!      a2=10.d0/1.d6
      b1=.45d0
!      b2=.45d0/1.d6
      b2=.45d0/1.d6*1.1d0
!      c1=0.03d0
	c1 = 0.00d0
!      c2=0.2d0/1.d12
      c2=0.1d0/1.d12
!      c3=4.5d6
      c3=4.0d6
      fphenom=a1*exp(-a2*qms)+b1*exp(-b2*qms)
!      fphenom=fphenom+c1*exp(-c2*(qms-4.5e6)**2)
      fphenom=fphenom+c1*exp(-c2*(qms-c3)**2)
      fphenom=sqrt(fphenom)
      return
      end

*fyukawa
      real*8 function fyukawa(qms,a)
      implicit real*8 (a-h,o-z)
      if(qms.lt.1.e-5.or.a.lt.1.e-5)then
      fyukawa=0.
      else
      arg=sqrt(qms/2.)/a
      fyukawa=atan(arg)/arg
      endif
      return
      end

*sigmot
      real*8 function sigmot(e,thr)
      implicit real*8 (a-h,o-z)
      alph=1./137.03604
      hbarc=197.3286
      sigmot=(alph*hbarc*cos(thr/2.)/2./e/sin(thr/2.)**2)**2
c  fm**2/sr
      return
      end

*recoil
      real*8 function recoil(e,thr,tm)
      implicit real*8 (a-h,o-z)
      recoil=1./(1.+2.*e*sin(thr/2.)**2/tm)
      return
      end

*sigqfs
      real*8 function sigqfs(e,th,w,z,a,eps,pf,delnu)
      implicit real*8 (a-h,o-z)
      character*1 dorg
      common /switch/dorg
c	Needed to save delnu for later calls to sigqfs
c	Saving of delsav between calls done by -fno-automatic

!      if(icall.ne.1) then
        delsav = delnu
!        icall = 1
!      else
!	delnu = delsav
!      end if
      pm=939.0d0
      up=2.7928456
      un=-1.91304184
      ap0=840.
!      ap1=750.
!      ap1=840.
      ap1=815.	! 5/2003
      alph=1./137.03604
      hbarc=197.32858
      pi=acos(-1.0d0)
      gamr=120.
      pfr=230.
      qmsrq=4.*730.*(730.-115.)*sin(37.1*pi/180./2.)**2
      qvsrq=qmsrq+115.**2
      na=int(a)
      if(na.eq.1)then
      ap=ap0
!      elseif(na.lt.4)then
!      ap=ap0+(a-1.)*(ap1-ap0)/3.
      elseif(na.lt.5)then
      ap=ap0+(a-1.)*(ap1-ap0)/4.
      else
      ap=ap1
      endif
c     print 200
  200 format(' enter de-e[MeV],domega-e[sr],b-luminosity[cm-2*s-1]')
c     read *,dee,dwe,blum
      thr=th*pi/180.0d0
      qms=4.0d0*e*(e-w)*sin(thr/2.0d0)**2
      qvs=qms+w**2
      ekappa=w-qms/2./pm
      if(ekappa.gt.-pm/2.)then
      cmtot=sqrt(pm**2+2.*pm*ekappa)
      else
      cmtot=pm
      endif
c  start qfs section
      signs=sigmot(e,thr)*recoil(e,thr,pm)
      call get_gegm_prot(qms,ap,gepsq,gmpsq)
      gmpsq=gmpsq*up**2
      formp=gepsq+qms*gmpsq/4./pm**2
      formp=formp/(1.+qms/4./pm**2)
      formp=formp+gmpsq*qms*tan(thr/2.)**2/2./pm**2
      sigep=signs*formp
      call get_gegm_neut(qms,ap,gensq,gmnsq)
      gmnsq=gmnsq*un**2
      formn=gensq+qms*gmnsq/4./pm**2
      formn=formn/(1.+qms/4./pm**2)
      formn=formn+gmnsq*qms*tan(thr/2.)**2/2./pm**2
      sigen=signs*formn
      epq=4.0d0*e**2*sin(thr/2.0d0)**2/2.0d0/pm
      epq=epq/(1.+2.*e*sin(thr/2.)**2/pm)+eps
      epq=e-epq
      if(int(a).eq.1)then

c	Get sig_el as ds2/(dE'dOmega) [nb/(sr MeV)]

	if(dorg.ne.'y') then	! Gauusian elastic peak

!      arg=(e-w-epq)/sqrt(2.)/1.	! Old qfs's wrong form of ds/dOmega
!      den=2.51
	  siggauss = delnu/2.354d0	! sigma_gauss = FWHM(=delnu)/2.354
	  showw = w
	  arg=(e-w-epq)/sqrt(2.0d0)/siggauss
	  den=2.51d0*siggauss		! 2.51 = sqrt(2*Pi)
	else

c	This is option for ds/dOmega

c	  wel=e - e/(1+2.0d0*e*sin(thr/2.)**2/pm)
c	  winv=pm*(pm+2.0d0*w)-qms
c 	  if(abs(w-wel).le.delnu.and.(w.le.wel))then
            sigq=(z*sigep+(a-z)*sigen)
c     	  else
c            sigq=0.0d0
c          endif
	  sigqfs=sigq
	  return
	end if	!	end dorg
      else
        gamq=gamr*pf*sqrt(qvs)/pfr/sqrt(qvsrq)
        arg=(e-w-epq)/1.20d0/(gamq/2.0d0)
        den=2.13d0*(gamq/2.d0)
!      sigq=(z*sigep+(a-z)*sigen)*exp(-arg**2)/den
      endif
      nq=int(arg)
!      if(abs(nq).gt.10)then
      if(abs(nq).gt.20)then	! needed for delw =1
        sigq=0.d0
      else
        sigq=(z*sigep+(a-z)*sigen)*exp(-arg**2)/den
      endif
      sigqfs=sigq
      return
      end

*sigdel
      real*8 function sigdel(e,th,w,a,epsd,pf,z)
      implicit real*8 (a-h,o-z)
      pm=939.0d0
      pimass=140.
      dm=1219.
!      ad1=700.
!      ad1=700.	! 2H
!      ad1=774.	! 12C
      ad1=750.	! 12C
!      ad0=700.
      ad0=774.
      pi=acos(-1.d0)
      alph=1./137.03604
      hbarc=197.32858
!      gamdp=110.
      gamdp=110.
!      gamsprd=140.
      gamsprd=20.
!      gamr=120.
      gamr=100.
!      gampi=5.
      gampi=50.
      qfdp=1.02d-7
      pfr=230.
      qmsr=4.*730.*(730.-390.)*sin(37.1*pi/180./2.)**2
      qvsr=qmsr+390.**2
      qmsrq=4.*730.*(730.-115.)*sin(37.1*pi/180./2.)**2
      qvsrq=qmsrq+115.**2
      na=int(a)
      if(na.eq.1)then
      qfd=qfdp
      gsprda=0.
      ad=ad0
!      elseif(na.lt.4)then
      elseif(na.eq.2)then
      qfd=qfdp
      gsprda=(a-1.)*gamsprd/3.
!      gsprda=gamspr
!      ad=ad0+(a-1.)*(ad1-ad0)/3.
      ad=700
      else
      ad=ad1
      gsprda=gamsprd
!      qfd=qfdp
      qfd=qfdp*1.1	! 5/03
      endif
      thr=th*pi/180.
      qms=4.*e*(e-w)*sin(thr/2.)**2
      qvs=qms+w**2
      ekappa=w-qms/2./pm
      cmtot2=pm**2+2.*pm*ekappa
c  begin delta calculation
      if(na.gt.1)then
      gamq=gamr*pf*sqrt(qvs)/pfr/sqrt(qvsrq)
      else
      gamq=0.
      endif
      epd=e-(dm-pm)*(dm+pm)/2./pm
      epd=epd/(1.+2.*e*sin(thr/2.)**2/pm)
      epd=epd-epsd
      wd=e-epd
      qmspk=4.*e*epd*sin(thr/2.)**2
      qvspk=qmspk+wd**2
c
c note width includes e-dependence,fermi broadening,& spreading
c
      wthresh=4.*e**2*sin(thr/2.)**2+pimass**2+2.*pimass*pm
      wthresh=wthresh/2./pm
       thrshfree=1.+2.*e*sin(thr/2.)**2/pm
      threshd=1.+pf/pm+pf**2/2./pm**2+2.*e*sin(thr/2.)**2/pm
       wthrfree=wthresh/thrshfree
      wthresh=wthresh/threshd
      if(w.gt.wthresh)then
        if((z.ne.1).and.(a.ne.1)) gampi = wthrfree-wthresh
      thresh=1.-exp(-(w-wthresh)/gampi)
      else
      thresh=0.
      endif
      gamd=gamdp
      gam=sqrt(gamd**2+gamq**2+gsprda**2)
!      sigd=qfdp*(gamdp/gam)
      sigd=qfd*(gamdp/gam)	! 5/03
      sigd=sigd*cmtot2*gam**2
      sigd=sigd/((cmtot2-(dm+epsd)**2)**2+cmtot2*gam**2)
      sigd=sigd*fd(qms,ad)**2/fd(qmsr,ad)**2
      test=qvs/qvsr
      sigd=sigd*test
      sigd=sigd*(qms/2./qvs+tan(thr/2.)**2)
      sigd=sigd/(qmsr/2./qvsr+tan(37.1*pi/180./2.)**2)
      sigd=sigd*sigmot(e,thr)/sigmot(730.d0,37.1*pi/180.)
      sigd=sigd*a                        
      sigd=sigd*thresh
      sigdel=sigd
      return
      end

*sigx
      real*8 function sigx(e,th,w,a)
      implicit real*8 (a-h,o-z)
      alph=1./137.03604
      pi=acos(-1.d0)
c     sig0=111.*1.e-30
      sig0=100.d-4
!      sig0=80.d-4
c     sig1=60.*1.e-27
      sig1=54.*1.d-1
      pimass=140.
      pm=939.0d0
c     gam0=550.
!     gam0=600.
      gam0=610.
!      gam0=700.
c     r=0.10
      aq=250.
      thr=th*pi/180.
      if(w.lt.1.e-5)go to 4
      qms=4.*e*(e-w)*sin(thr/2.)**2
      arg0=w-qms/2./pm-pimass-pimass**2/2./pm
      arg1=arg0/gam0
      arg=arg1**2/2.
      if(arg1.gt.8.)then
      shape=1.+sig1/sig0/arg0
      elseif(arg1.lt.1.e-5)then
      shape=0.
      elseif(arg1.lt.0.1)then
      shape=sig1*arg0/2./gam0**2/sig0
      else
      shape=(1.-exp(-arg))*(1.+sig1/sig0/arg0)
      endif
      ekappa=w-qms/2./pm
      siggam=sig0*shape
      qs=qms+w**2
      eps=1./(1.+2.*qs*tan(thr/2.)**2/qms)
      flux=alph*ekappa*(e-w)/2./pi**2/qms/e/(1.-eps)
      if(flux.lt.1.e-20)flux=0.
      sigee=flux*siggam*fphenom(qms)**2
c     sigee=flux*siggam
!      r=0.56*1.e6/(qms+pm**2)
	r=min(0.56*1.e6/(qms+pm**2),0.2)
      factor1=1.+eps*r
      sigee=sigee*factor1
 4    sigx=a*sigee
      return
      end

*sigr1
      real*8 function sigr1(e,th,w,a,pf,z)
      implicit real*8 (a-h,o-z)
      pi=acos(-1.d0)
      pm=939.0d0
      pimass=140.
      thr=th*pi/180.
      pfr=230.
!     rm=1500.
      rm=1500.
      epsr=0.
      ar0=1000.
!      ar1=1000.
!      ar1=940.	! 2H
      ar1=1000.	! 12C
!      gamqfr=120.
      gamqfr=100.
!      gamsprd=140.
      gamsprd=0.
      gamr=110.
!      gamr=130.
!      gampi=5.
      gampi=25.
      qfrp=1.20d-7
      qmsqfr=4.*730.*(730.-115.)*sin(37.1*pi/180./2.)**2
      qvsqfr=qmsqfr+115.**2
      qmsrr=4.*10000.*(10000.-1240.)*sin(6.*pi/180./2.)**2
      qvsrr=qmsrr+1240.**2
      sigref=fd(qmsrr,ar0)**2*qvsrr
      sigref=sigref*(qmsrr/2./qvsrr+tan(6.*pi/180./2.)**2)
      sigref=sigref*sigmot(10000.d0,6.*pi/180.)
      na=int(a)
      if(na.eq.1)then
      qfr=qfrp
      gsprda=0.
      ar=ar0
!      elseif(na.lt.4)then
      elseif(na.eq.2)then
      qfr=qfrp
      gsprda=(a-1.)*gamsprd/3.
!      ar=ar0+(a-1.)*(ar1-ar0)/3.
      ar=920
      else
      ar=ar1
      gsprda=gamsprd
      qfr=qfrp
      endif
      qms=4.*e*(e-w)*sin(thr/2.)**2
      qvs=qms+w**2
      if(na.gt.1)then
      gamq=gamqfr*pf*sqrt(qvs)/pfr/sqrt(qvsqfr)
      else
      gamq=0.
      endif
      cmtot2=pm**2+2.*pm*w-qms
      wthresh=4.*e**2*sin(thr/2.)**2+pimass**2+2.*pimass*pm
      wthresh=wthresh/2./pm
       thrshfree=1.+2.*e*sin(thr/2.)**2/pm
      threshd=1.+pf/pm+pf**2/2./pm**2+2.*e*sin(thr/2.)**2/pm
       wthrfree=wthresh/thrshfree
      wthresh=wthresh/threshd
      if(w.gt.wthresh)then
        if((z.ne.1).and.(a.ne.1)) gampi = wthrfree-wthresh
      thresh=1.-exp(-(w-wthresh)/gampi)
      else
      thresh=0.
      endif
      epr=e-(rm-pm)*(rm+pm)/2./pm
      epr=epr/(1.+2.*e*sin(thr/2.)**2/pm)
      epr=epr-epsr
      wr=e-epr
      gam=sqrt(gamr**2+gamq**2+gsprda**2)
      sigr=qfr*(gamr/gam)/sigref
      sigr=sigr*cmtot2*gam**2
      sigr=sigr/((cmtot2-(rm+epsr)**2)**2+cmtot2*gam**2)
      sigr=sigr*qvs*fd(qms,ar)**2
      sigr=sigr*(qms/2./qvs+tan(thr/2.)**2)
      sigr=sigr*sigmot(e,thr)
      sigr1=a*thresh*sigr
      return
      end

*sigr2
      real*8 function sigr2(e,th,w,a,pf,z)
      implicit real*8 (a-h,o-z)
      pi=acos(-1.d0)
      pm=939.0d0
      pimass=140.
      thr=th*pi/180.
      pfr=230.
      rm=1700.
      epsr=0.
      ar0=1200.
!      ar1=1200.
!      ar1=1000.	! 2H
      ar1=1200.	! 12C
!      gamqfr=120.
      gamqfr=100.
!      gamsprd=140.
      gamsprd=0.
      gamr=110.
!      gampi=5.
      gampi=25.
      qfrp=0.68d-7
      qmsqfr=4.*730.*(730.-115.)*sin(37.1*pi/180./2.)**2
      qvsqfr=qmsqfr+115.**2
      qmsrr=4.*10000.*(10000.-1520.)*sin(6.*pi/180./2.)**2
      qvsrr=qmsrr+1520.**2
      sigref=fd(qmsrr,ar0)**2*qvsrr
      sigref=sigref*(qmsrr/2./qvsrr+tan(6.*pi/180./2.)**2)
      sigref=sigref*sigmot(10000.d0,6.*pi/180.)
      na=int(a)
      if(na.eq.1)then
      qfr=qfrp
      gsprda=0.
      ar=ar0
!      elseif(na.lt.4)then
      elseif(na.eq.2)then
      qfr=qfrp
      gsprda=(a-1.)*gamsprd/3.
!      ar=ar0+(a-1.)*(ar1-ar0)/3.
	 ar=980
      else
      ar=ar1
      gsprda=gamsprd
      qfr=qfrp
      endif
      qms=4.*e*(e-w)*sin(thr/2.)**2
      qvs=qms+w**2
      if(na.gt.1)then
      gamq=gamqfr*pf*sqrt(qvs)/pfr/sqrt(qvsqfr)
      else
      gamq=0.
      endif
      cmtot2=pm**2+2.*pm*w-qms
      wthresh=4.*e**2*sin(thr/2.)**2+pimass**2+2.*pimass*pm
      wthresh=wthresh/2./pm
       thrshfree=1.+2.*e*sin(thr/2.)**2/pm
      threshd=1.+pf/pm+pf**2/2./pm**2+2.*e*sin(thr/2.)**2/pm
       wthrfree=wthresh/thrshfree
      wthresh=wthresh/threshd
      if(w.gt.wthresh)then
        if((z.ne.1).and.(a.ne.1)) gampi = wthrfree-wthresh
      thresh=1.-exp(-(w-wthresh)/gampi)
      else
      thresh=0.
      endif
      epr=e-(rm-pm)*(rm+pm)/2./pm
      epr=epr/(1.+2.*e*sin(thr/2.)**2/pm)
      epr=epr-epsr
      wr=e-epr
      gam=sqrt(gamr**2+gamq**2+gsprda**2)
      sigr=qfr*(gamr/gam)/sigref
      sigr=sigr*cmtot2*gam**2
      sigr=sigr/((cmtot2-(rm+epsr)**2)**2+cmtot2*gam**2)
      sigr=sigr*qvs*fd(qms,ar)**2
      sigr=sigr*(qms/2./qvs+tan(thr/2.)**2)
      sigr=sigr*sigmot(e,thr)
      sigr2=a*thresh*sigr
      return
      end

*sig2n
      real*8 function sig2n(e,th,w,z,a,pf)
      implicit real*8 (a-h,o-z)
      pi=acos(-1.d0)
      thr=th*pi/180.d0
      dm=1219.d0
      pimass=140.d0
      pm=939.d0
!      a2=550.
      a2=500.d0
      pfr=60.d0
      gam2n=20.d0
!      gam2n=40.
      gamqfr=40.d0
      gamref=300.d0
      gamr=gamref
      sigref=0.20d-7
      qmsr=4.d0*596.8d0*(596.8d0-380.d0)*sin(60.d0*pi/180.d0/2.d0)**2
      qvsr=qmsr+380.d0**2
      sigkin=0.5d0*sigmot(596.8d0,60.d0*pi/180.d0)
      sigkin=sigkin*(qmsr/2.d0/qvsr+tan(60.d0*pi/180.d0/2.d0)**2)
      sigkin=sigkin*qvsr*fd(qmsr,a2)**2
      sigkin=sigkin*gamr/gamref
      sigcon=sigref/sigkin
      qms=4.d0*e*(e-w)*sin(thr/2.d0)**2
      qvs=qms+w**2
      gamqf=gamqfr*(pf/pfr)*(sqrt(qvs)/sqrt(qvsr))
      effmass=(pm+dm)/2.d0
      sig=(z*(a-z)/a)*sigmot(e,thr)
      sig=sig*(qms/2.d0/qvs+tan(thr/2.d0)**2)
      sig=sig*qvs*fd(qms,a2)**2
      ekappa=w-qms/2.d0/pm
      cmtot2=pm**2+2.d0*pm*ekappa
c     gam=sqrt(gamr**2+gamqf**2)
      gam=gamr
      sig=sig*cmtot2*gam**2
      sig=sig/((cmtot2-effmass**2)**2+cmtot2*gam**2)
      sig=sig*(gamr/gam)*sigcon
      sig2n=sig
      wthresh=qms/4.d0/pm
      if(w.gt.wthresh)then
      thresh=1.d0-exp(-(w-wthresh)/gam2n)
      else
      thresh=0.d0
      endif
      sig2n=sig2n*thresh
      return
      end

*radiate
      subroutine radiate(e,th,w,z,a,signr,sigrad)
      implicit real*8 (a-h,o-z)
      common /par/ eps,epsd,pf,spence
      alph=1./137.03604
      emass=0.511
      del=10.
      pi=acos(-1.d0)
      prec=.0005
      thr=th*pi/180.
      arg=cos(thr/2.)**2
      spence=pi**2/6.-log(arg)*log(1.-arg)
      do 10 nsp=1,50
 10   spence=spence-arg**nsp/float(nsp)**2
c     print 15,spence
 15   format(' spence function = ',1f14.9)
      qms=4.*e*(e-w)*sin(thr/2.)**2
      d1=(2.*alph/pi)*(log(qms/emass**2)-1.)
      d2=13.*(log(qms/emass**2)-1.)/12.-17./36.-0.5*(pi**2/6.-spence)
      d2=d2*(2.*alph/pi)
      ebar=sqrt(e*(e-w))
      sigrad=signr*(1.+d2)*exp(-d1*log(ebar/del))
c     print 20,d1,d2
 20   format(' delta1 and delta2 = ',2e14.6)
      x1=0.
      x2=w-del
      call rom(e,th,w,z,a,x1,x2,prec,ans,kf)
      ans=ans*1.d-26
      iflag=0
      err=0.
c     print 25,x1,x2,prec,ans,err
 25   format(' x1,x2,prec,ans,err = ',5e10.3)
c     print 30,kf,iflag
 30   format(' kf,iflag = ',2i5)
      sigrad=sigrad+ans
      return
      end

*rom
      subroutine rom(e0,th0,w0,z0,a0,a,b,eps,ans,k)
      implicit real*8 (a-h,o-z)
c  romberg method of integration
      dimension w(50,50)
      h=b-a
      k=0
      call value(e0,th0,w0,z0,a0,a,fa)
      call value(e0,th0,w0,z0,a0,b,fb)
      w(1,1)=(fa+fb)*h/2.
    4 k=k+1
      if(k.ge.49)go to 5
      h=h/2.
      sig=0.
      m=2**(k-1)
      do 1 j=1,m
      j1=2*j-1
      x=a+float(j1)*h
      call value(e0,th0,w0,z0,a0,x,f)
    1 sig=sig+f
      w(k+1,1)=w(k,1)/2.+h*sig
      do 2 l=1,k
      iu=k+1-l
      iv=l+1
    2 w(iu,iv)=(4.**(iv-1)*w(iu+1,iv-1)-w(iu,iv-1))/(4.**(iv-1)-1.)
      e=(w(iu,iv)-w(iu,iv-1))/w(iu,iv)
      if(abs(e)-eps) 3,3,4
    3 ans=w(1,iv)
      return
    5 print 100
  100 format(' k overflow')
      stop
!      call exit
      end

*value
      subroutine value(e,th,w,z,a,x,f)
      implicit real*8 (a-h,o-z)
      common /par/ eps,epsd,pf,spence
      alph=1./137.03604
      emass=.511
      pi=acos(-1.d0)
      thr=th*pi/180.
      sig1=sigqfs(e,th,x,z,a,eps,pf,delnu)
      sig1=sig1+sigdel(e,th,x,a,epsd,pf,z)
      sig1=sig1+sigx(e,th,x,a)
      sig1=sig1+sigr1(e,th,x,a,pf,z)
      sig1=sig1+sigr2(e,th,x,a,pf,z)
      sig1=sig1+sig2n(e,th,x,z,a,pf)
      sig2=sigqfs(e-w+x,th,x,z,a,eps,pf,delnu)
      sig2=sig2+sigdel(e-w+x,th,x,a,epsd,pf,z)
      sig2=sig2+sigx(e-w+x,th,x,a)
      sig2=sig2+sigr1(e-w+x,th,x,a,pf,z)
      sig2=sig2+sigr2(e-w+x,th,x,a,pf,z)
      sig2=sig2+sig2n(e-w+x,th,x,z,a,pf)
      qms1=4.*e*(e-x)*sin(thr/2.)**2
      qms2=4.*(e-w+x)*(e-w)*sin(thr/2.)**2
      qmsbar=sqrt(qms1*qms2)
      f1=(log(qms1/emass**2)-1.)/2./(e-x)**2
      f1=f1*sig1
      f1=f1*alph*((e-x)**2+(e-w)**2)/(w-x)/pi
      f2=sig2*(log(qms2/emass**2)-1.)/2./e**2
      f2=f2*alph*(e**2+(e-w+x)**2)/(w-x)/pi
      d1=(2.*alph/pi)*(log(qmsbar/emass**2)-1.)
      d2=13.*(log(qmsbar/emass**2)-1.)/12.-17./36.
      d2=(2.*alph/pi)*(d2-0.5*(pi**2/6.-spence))
      ebar=sqrt((e-x)*(e-w))
      f=(f1+f2)*(1.+d2)*((w-x)/ebar)**d1
      return
      end
*
      subroutine get_gegm_prot(q2,ap,gepsq,gmpsq)
      real*8 q2,ap,gepsq,gmpsq,fd,gm,ge
      character*8 model
      real*8 fd1,fda,dipcorr
c
c   return Ge^2 and Gm^2/u  
c
      model = 'bosted'
      if ( model .eq. 'dipole') then
      gepsq = fd(q2,ap)**2
      gmpsq = fd(q2,ap)**2
      endif
      if ( model .eq. 'bosted') then

c A>1 dipole correction (bound nucleons)
	fd1 = fd(q2,840.d0)
	fda = fd(q2,ap)
	dipcorr = fda/fd1

         q2 = q2/1.d6
        gm=1.d0+0.35d0*sqrt(q2)+2.44d0*q2+.5d0*sqrt(q2)*q2
         gm = 1.d0/(gm +1.04d0*q2*q2 + .34d0*q2*q2*sqrt(q2))
         ge=1.d0/(1.d0+0.62d0*sqrt(q2)+0.68d0*q2
	1	+2.8d0*sqrt(q2)*q2+0.84d0*q2*q2)

c A>1 dipole correction (bound nucleons)
	ge = ge*dipcorr
	gm = gm*dipcorr
	
      gepsq = ge**2
      gmpsq = gm**2
      q2 = 1.d6*q2
      endif
!      write(*,*) q2/1000./1000.,gepsq,gmpsq
      return
      end
c
c
c
c
      subroutine get_gegm_neut(q2,ap,gensq,gmnsq)
      real*8 q2,ap,gensq,gmnsq,fd,gm,ge,gb,mp,falloff,un
      character*8 model
      real*8 fd1,fda,dipcorr
c
c   return Ge^2 and Gm^2/u  
c
c model = 'walcher' is from J. Friedrich and Th. Walcher fit
c               in Eur.Phys.J. A17 (2003) 607-623
      model = 'walcher'
      if ( model .eq. 'dipole') then
      mp=939.d0 ! as done in sigqfs
      falloff=(1.d0+5.6d0*q2/4.d0/mp**2)**2
      un=-1.91304184d0
      gensq = fd(q2,ap)**2/falloff*(un*q2/4.d0/mp**2)**2
      gmnsq = fd(q2,ap)**2
      endif
      if ( model .eq. 'walcher') then

c A>1 dipole correction (bound nucleons)
	fd1 = fd(q2,840.d0)
	fda = fd(q2,ap)
	dipcorr = fda/fd1

         q2 = q2/1.d6
      ge = 1.04d0/(1.d0+q2/1.73d0)**2 - 1.04d0/(1.d0+q2/1.54d0)**2
      gb = exp(-1.d0/2.d0*((sqrt(q2)-.29d0)/.20d0)**2) 
     >     +  exp(-1.d0/2.d0*((sqrt(q2)+.29d0)/.20d0)**2)
      ge = ge + 0.23d0*q2*gb
      gm = 1.012d0/(1.d0+q2/0.77d0)**2 - 0.012d0/(1+q2/6.8d0)**2
      gb = exp(-1./2.*((sqrt(q2)-.33)/.14)**2) 
     >     +  exp(-1./2.*((sqrt(q2)+.33)/.14)**2)
      gm = gm - 0.28*q2*gb
      gd = 1/(1+q2/.71)/(1+q2/.71)
      q2 = 1000.*1000.*q2

c A>1 dipole correction (bound nucleons)
	ge = ge*dipcorr
	gm = gm*dipcorr

      gensq = ge**2
      gmnsq = gm**2
      endif
      return
      end
