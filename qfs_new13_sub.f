*qfs
!      subroutine qfs(z,a,e,th,w,delnu,flag,sigma_qfs,elflag)
      subroutine qfs(z,a,e,th,w,delnu,flag,sigma_qfs,elflag,scalein)! OR 3/05
cc    subroutine qfs(z,a,e,th,w,delnu,flag,elflag,sigma_qfs) 
c
      implicit none
      real*8 z,a,e,th,w,delnu,sigma_qfs,elflag,scalein
      real*8 scale, wpi, eps,epsd,pf,spence, pm, dm, alph, pi, sigqfza,
     &       sigqfs, thr, qsq, wsq, f1, f2, r, siginel, sigmot, sigdelc,
     &       sigr1c, sigr2c, sigxc, sigda, sigr1a, sigr2a, sigxa, 
     &       sig2na, sigdel, sigx, sigr1, sigr2, sig2n, sig, sigrad,
     &       escat,xsecqe,xsecdis,xsectot

ccc   Fro convenience, define and use these vars (used for deut inel xsec)
      real*8 sig_mott, alpha, W1, W2, Mp, e0, ep, e_nu, hbarcsq, sigma_qe
      real*8 scale_peterxsec
      real*8 r09
 
      character*1 dorg
      logical*4 flag
      character*11 proton_inel_fit

      common /par/ eps,epsd,pf,spence
      common /switch/ dorg
      parameter(scale_peterxsec=1.d-36)  !! multiply this to peter's xsec
      parameter(scale=1.d-26)
      parameter (wpi = 1.078d0)
	
      pm=939.0d0
      dm=1219.d0
      alph=1.d0/137.03604d0
      pi=acos(-1.d0)
c
      alpha = 0.007297352533d0 !! define or reassign some vars for convenience
      Mp  = pm/ 1000.0d0       !!  these are in GeV !!! 
      e0  = e / 1000.0d0               
      e_nu= w / 1000.0d0
      ep = e0 - e_nu
      hbarcsq=0.389379292d0
c     print *, 'e_nu=', e_nu

c
c   Select which Hall C proton inelastic fit to use.
c 
c        proton_inel_fit = 'liang      '
c       proton_inel_fit = 'EC_nov_2004'
c        proton_inel_fit = 'EC_050614'
        proton_inel_fit = 'PB_MEC_09'
c
      IF(A.EQ.1.d0) THEN
        PF=0.d0    
        EPS=0.d0
        EPSD=13.0d0
        dorg='y'
      ELSEIF(A.eq.2.d0) THEN
!        PF=77.d0 
!        EPS=2.23d0
!        PF=74.d0	!	OR 2/04
        PF=100.d0	!	OR 2/05
        EPS=2.22d0
        EPSD=13.0d0
      ELSEIF(A.EQ.4.d0) THEN
!        PF=180.d0  !152   
!        EPS=20.d0  !8
        PF=180.d0  ! 	OR 3/05   
        EPS=20.2d0  !	OR 3/05
cc      PF=162.d0  ! 	OR 3/05   
cc      EPS=20.6d0  !	OR 3/05
        EPSD=13.0d0
      ELSEIF(A.EQ.9.d0) THEN
        PF=200.d0       
        EPS=20.d0
        EPSD=13.0d0
      ELSEIF(A.EQ.12.d0) THEN
        PF=230.d0    !205, 215, 220, 230=5/03     
        EPS=25.d0 
        EPSD=13.d0
      ELSEIF(A.EQ.14.d0 .OR. A.EQ.15.d0) THEN
        PF=225.d0      
        EPS=25.d0
        EPSD=13.d0
      ELSEIF(A.EQ.27.d0) THEN
        PF=235.d0       
        EPS=32.d0       
        EPSD=13.d0
      ELSEIF(A.EQ.59.d0) THEN
        PF=260.d0               
        EPS=35.d0
        EPSD=13.d0
      ELSEIF(A.EQ.64.d0) THEN
        PF=260.d0                
        EPS=35.d0
        EPSD=13.d0
      ELSEIF(A.EQ.184.d0) THEN
        PF=265.d0                
        EPS=43.d0
        EPSD=13.d0
      ELSE
        WRITE(*,*) 'pf, eps, epsd not specified for this material, stop'
        STOP
      ENDIF
c
!	6/03 Use Hall C F2 R for proton inelastic sigma

      thr = th*pi/180.0d0	! 6/03
      qsq = 4.0d0*e*(e-w)*sin(thr/2.0d0)**2	! 6/03
      wsq = pm*pm + 2.d0*pm*w - qsq	! 6/03

      if((a.eq.1.d0 .and.z.eq.1.d0)) then
         siginel = 0.
         sigqfza =0.
         if ( elflag.eq.0.d0 .and.(wsq/1.d6.le.9.0d0)) then	! 6/03

          if ( proton_inel_fit .eq. 'liang      ') then
	   call hcf2r(qsq/1.d6,wsq/1.d6,f2,r)	! 6/03
          else if ( proton_inel_fit .eq. 'EC_nov_2004') then
           call christy_rss(qsq/1.d6,wsq/1.d6,f2,r)
          else if ( proton_inel_fit .eq. 'EC_050614') then
           call ressf(qsq/1.d6,wsq/1.d6,f2,r)
          else if ( proton_inel_fit .eq. 'PB_MEC_09') then
           call F1F2IN09(1.0d0, 1.0d0,qsq/1.d6,wsq/1.d6, F1, F2, r)
          endif
	   siginel=sigmot(e,thr)*f2*(1.d0/w+
     &     2.d0*w*(1.d0+qsq/w**2)*(tan(thr/2.d0))**2/(qsq*(1d0+r))) ! 6/03
	   siginel=siginel*scale	! 6/03
c           write(*,*) qsq/1.d6,wsq/1.d6,siginel
         else
          sigqfza=sigqfs(e,th,w,z,a,eps,pf,delnu)*scale
         endif
      else				! 6/03 A >1

        if(a.eq.2.d0 .and.z.eq.1.d0) then	! 2/04

!	  print *,'e ',e,enu,ww,qsq,th
c		call hallc2h(e,th,w,sigdelc,sigr1c,sigr2c,sigxc)	! 2/04
c		sigda = sigdelc*scale	!	"
c		sigr1a = sigr1c*scale	!	"
c		sigr2a = sigr2c*scale	!	"
c		sigxa = sigxc*scale		!	"
c	       siginel = sigda+sigxa+sigr1a+sigr2a	!	"
c

c              print *, 'before call f1f2in06:', z, a, qsq/1.d6, wsq/1.d6
c               call F1F2IN06(Z, A, qsq/1.d6, wsq/1.d6, F1, F2)  !!peter's code
               call F1F2IN09(Z, A, qsq/1.d6, wsq/1.d6, F1, F2,R09)  !! Peter and Vahe new code 09
c              print *, 'after  call f1f2in06:', F1, F2
               W1 = F1 / Mp
               W2 = F2 / e_nu
               sig_mott = alpha**2 / (qsq/1.d6) / tan(thr/2.0d0)**2 * ep/e0
               sig_mott = sig_mott * hbarcsq * 1.0d6  !!xsec in nb/GeV/str/nuc
               siginel = sig_mott * ( W2 + 2.0d0*W1*tan(thr/2.0d0)**2 )
               siginel = siginel * scale_peterxsec
c              print *,'Winv=',sqrt(wsq),' Q2=',(qsq/1.d6),' siginel=',siginel

		sig2na = 0.d0	!	"

cccccccccccc   ***** y-scaling routine ******** cccccccc
cc 
cc     add call to y-scaling fit to e-D quasi-free to get quasi-free xn
cc      qfs_deut expects energy in GeV, remember that "w" is nu=e-escat 
cc          returns the quasi-elastic xn, inelastic xn and the sum.
cc          Will only use the quasi-elastic xn
c                escat = (e-w)/1000.d0
c                call qfs_deut(e/1000.d0,thr,escat,xsecqe,xsecdis,xsectot)
cc                write(*,*) ' comp q.e. xn ',sigqfza,xsecqe*scale
c             sigqfza = xsecqe*scale ! comment out this line to use old version
ccccccccccc

cc         Let's use peter's QE routine: 

c                call QUASIY8_2006(Z,A,E0,EP,TH,SIGMA_qe)  !! "_2006" added
c new for 09
                call F1F2QE09(Z, A, qsq/1.d6, wsq/1.d6, F1, F2)
               W1 = F1 / Mp
               W2 = F2 / e_nu
               sig_mott = alpha**2 / (qsq/1.d6) / tan(thr/2.0d0)**2 * ep/e0
               sig_mott = sig_mott * hbarcsq * 1.0d6  !!xsec in nb/GeV/str/nuc
               sigma_qe = sig_mott * ( W2 + 2.0d0*W1*tan(thr/2.0d0)**2 )
c end new for 09
                sigma_qe = sigma_qe * scale_peterxsec 
                sigqfza = sigma_qe 
c               print *, "sigma_qe=", sigma_qe

	 else	!	A > 2

         sigda=sigdel(e,th,w,a,epsd,pf,z)*scale
         sigxa=sigx(e,th,w,a)*scale
         sigr1a=sigr1(e,th,w,a,pf,z)*scale
         sigr2a=sigr2(e,th,w,a,pf,z)*scale
!         siginel=sigda+sigxa+sigr1a+sigr2a+sig2na	! 6/03
	  siginel=sigda+sigxa+sigr1a+sigr2a

   	  siginel=siginel*scalein	! OR 3/05: scalein from *.inp file, inelastic scale

!	  sig2na=sig2n(e,th,enu,z,a,pf)*scale*0.7d0
         sig2na=sig2n(e,th,w,z,a,pf)*scale*1.0d0		! 11/04
	  siginel=siginel+sig2na
c add calls to Peter's code for A> 2 mkj Apr 16 2007
c               call F1F2IN06(Z, A, qsq/1.d6, wsq/1.d6, F1, F2)  !!peter's code
               call F1F2IN09(Z, A, qsq/1.d6, wsq/1.d6, F1, F2,R09)  !! Peter and Vahe new code 09
c              print *, 'after  call f1f2in06:', F1, F2
               W1 = F1 / Mp
               W2 = F2 / e_nu
               sig_mott = alpha**2 / (qsq/1.d6) / tan(thr/2.0d0)**2 * ep/e0
               sig_mott = sig_mott * hbarcsq * 1.0d6  !!xsec in nb/GeV/str/nuc
               siginel = sig_mott * ( W2 + 2.0d0*W1*tan(thr/2.0d0)**2 )
               siginel = siginel * scale_peterxsec*scalein 
c              call QUASIY8_2006(Z,A,E0,EP,TH,SIGMA_qe)  !! "_2006" added
c new for 09
                call F1F2QE09(Z, A, qsq/1.d6, wsq/1.d6, F1, F2)
               W1 = F1 / Mp
               W2 = F2 / e_nu
               sig_mott = alpha**2 / (qsq/1.d6) / tan(thr/2.0d0)**2 * ep/e0
               sig_mott = sig_mott * hbarcsq * 1.0d6  !!xsec in nb/GeV/str/nuc
               sigma_qe = sig_mott * ( W2 + 2.0d0*W1*tan(thr/2.0d0)**2 )
c end new for 09
               sigma_qe = sigma_qe * scale_peterxsec 
                sigqfza = sigma_qe 
        
c
	 end if ! hallc 2h

!	   sig2na=sig2n(e,th,enu,z,a,pf)*scale*0.7d0
c  
 
	end if	! hcf2r	2/04
!      end if						! 6/03
      
!      sig=sigqfza+sigda+sigxa+sigr1a+sigr2a+sig2na
	sig = sigqfza + siginel				! 6/03
c
        if ( (a.eq.1.d0 .and. z.eq.1.d0) ) then
         if (  elflag .eq. 1.d0) then 
          sig = sigqfza
         else
           sig = siginel 
         endif
        endif
c	  
      sigma_qfs=sig        ! nonradiated qfs cross section
c
c      thcr=th*pi/180.d0
c      qms=4.d0*e*(e-w)*sin(thr/2.d0)**2
c      qvs=qms+w**2
c      ekappa=w-qms/2.d0/pm
c      if(ekappa.gt.-pm/2.d0)then
c        cmtot=sqrt(pm**2+2.d0*pm*ekappa)
c      else
c        cmtot=pm
c      endif
c      flux=(alph/2.d0/pi**2)*((e-w)/e)*((2.d0*pm*w-qms)/2.d0/pm/qms)
c      polari=1.d0/(1.d0+2.d0*qvs*tan(thr/2.d0)**2/qms)
c      flux=flux/(1.d0-polari)
c      if(ekappa.lt.0.) photsig=0.d0
c      photsig=sig/flux
c radiate
      sigrad=0.d0
      if(sig.gt.0.0d0.and.flag) then
	call radiate(e,th,w,z,a,sig,sigrad,delnu)  ! 'delnu' added mar05
        sigma_qfs=sigrad
      end if
c
      return
      end

*fd
      real*8 function fd(qms,a)
      implicit none
      real*8 qms, a
      fd=1.d0/(1.d0+qms/a**2)**2
      return
      end

*fm
      real*8 function fm(qms,a)
      implicit none
      real*8 qms, a
      fm=1.d0/(1.d0+qms/a**2)
      return
      end

*fphenom
      real*8 function fphenom(qms)
      implicit none
      real*8 qms, a1, a2, b1, b2, c1, c2, c3
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
      fphenom=a1*exp((-a2)*qms)+b1*exp((-b2)*qms)
!      fphenom=fphenom+c1*exp((-c2)*(qms-4.5d6)**2)
      fphenom=fphenom+c1*exp((-c2)*(qms-c3)**2)
      fphenom=sqrt(fphenom)
      return
      end

*fyukawa
      real*8 function fyukawa(qms,a)
      implicit none
      real*8 qms, a, arg
      if(qms.lt.1.d-5.or.a.lt.1.d-5)then
      fyukawa=0.d0
      else
      arg=sqrt(qms/2.d0)/a
      fyukawa=atan(arg)/arg
      endif
      return
      end

*sigmot
      real*8 function sigmot(e,thr)
      implicit none
      real*8 e, thr, alph, hbarc
      alph=1.d0/137.03604d0
      hbarc=197.3286d0
      sigmot=(alph*hbarc*cos(thr/2.d0)/2.d0/e/sin(thr/2.d0)**2)**2
c  fm**2/sr
      return
      end

*recoil
      real*8 function recoil(e,thr,tm)
      implicit none
      real*8 e,thr,tm
      recoil=1.d0/(1.d0+2.d0*e*sin(thr/2.d0)**2/tm)
      return
      end

*sigqfs
      real*8 function sigqfs(e,th,w,z,a,eps,pf,delnu)
      implicit none
      real*8 e,th,w,z,a,eps,pf,delnu
      real*8  pm, up, un, ap0, ap1, alph, hbarc, pi, gamr, pfr,
     &       qmsrq, qvsrq, ap, thr, qms, qvs, ekappa, cmtot, signs, 
     &       sigmot, recoil, gepsq, gmpsq, formp, sigep, gensq, gmnsq,
     &       formn, sigen, epq, siggauss, den, arg, sigq, gamq
      integer na, nq
 
      character*1 dorg
      common /switch/dorg
c	Needed to save delnu for later calls to sigqfs
c	Saving of delsav between calls done by -fno-automatic

!      if(icall.ne.1) then
!        delsav = delnu
!        icall = 1
!      else
!	delnu = delsav
!      end if
      pm=939.0d0
      up=2.7928456d0
      un=-1.91304184d0
      ap0=840.d0
!      ap1=750.d0
!      ap1=840.d0
      ap1=815.d0	! 5/2003
      alph=1.d0/137.03604d0
      hbarc=197.32858d0
      pi=acos(-1.0d0)
      gamr=120.d0
      pfr=230.d0
      qmsrq=4.d0*730.d0*(730.d0-115.d0)*sin(37.1d0*pi/180.d0/2.0d0)**2
      qvsrq=qmsrq+115.d0**2
      na=int(a)
      if(na.eq.1)then
      ap=ap0
!      elseif(na.lt.4)then
!      ap=ap0+(a-1.d0)*(ap1-ap0)/3.d0
      elseif(na.lt.5)then
      ap=ap0+(a-1.d0)*(ap1-ap0)/4.d0
      else
      ap=ap1
      endif
c     print 200
  200 format(' enter de-e[MeV],domega-e[sr],b-luminosity[cm-2*s-1]')
c     read *,dee,dwe,blum
      thr=th*pi/180.0d0
      qms=4.0d0*e*(e-w)*sin(thr/2.0d0)**2
      qvs=qms+w**2
      ekappa=w-qms/2.d0/pm
      if(ekappa.gt.(-pm)/2.d0)then
      cmtot=sqrt(pm**2+2.d0*pm*ekappa)
      else
      cmtot=pm
      endif
c  start qfs section
      signs=sigmot(e,thr)*recoil(e,thr,pm)
      call get_gegm_prot(qms,ap,gepsq,gmpsq)
      gmpsq=gmpsq*up**2
      formp=gepsq+qms*gmpsq/4.d0/pm**2
      formp=formp/(1.d0+qms/4.d0/pm**2)
      formp=formp+gmpsq*qms*tan(thr/2.d0)**2/2.d0/pm**2
      sigep=signs*formp
      call get_gegm_neut(qms,ap,gensq,gmnsq)
      gmnsq=gmnsq*un**2
      formn=gensq+qms*gmnsq/4.d0/pm**2
      formn=formn/(1.d0+qms/4.d0/pm**2)
      formn=formn+gmnsq*qms*tan(thr/2.d0)**2/2.d0/pm**2
      sigen=signs*formn
      epq=4.0d0*e**2*sin(thr/2.0d0)**2/2.0d0/pm
      epq=epq/(1.d0+2.d0*e*sin(thr/2.d0)**2/pm)+eps
      epq=e-epq
      if(int(a).eq.1)then

c	Get sig_el as ds2/(dE'dOmega) [nb/(sr MeV)]

	if(dorg.ne.'y') then	! Gauusian elastic peak

!      arg=(e-w-epq)/sqrt(2.)/1.d0	! Old qfs's wrong form of ds/dOmega
!      den=2.51d0
	  siggauss = delnu/2.354d0	! sigma_gauss = FWHM(=delnu)/2.354
c	  showw = w
	  arg=(e-w-epq)/sqrt(2.0d0)/siggauss
	  den=2.51d0*siggauss		! 2.51 = sqrt(2*Pi)
	else

c	This is option for ds/dOmega

c	  wel=e - e/(1.d0+2.0d0*e*sin(thr/2.d0)**2/pm)
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
        sigq=(z*sigep+(a-z)*sigen)*exp(-(arg**2))/den
      endif
      sigqfs=sigq
      return
      end

*sigdel
      real*8 function sigdel(e,th,w,a,epsd,pf,z)
      implicit none
      real*8 e,th,w,a,epsd,pf,z
      real*8 alph, pi, pimass, pm, thr, qms, ekappa, dm, ad1, ad0, 
     &       hbarc, gamdp, gamsprd, gamr, gampi, qfdp, pfr, qmsr,
     &       qvsr, qmsrq, qvsrq, qfd, gsprda, ad, qvs, cmtot2,
     &       gamq, epd, wd, qmspk, wthresh, thrshfree, threshd,
     &       wthrfree, thresh, gamd, gam, sigd, fd, test, 
     &       sigmot
      integer na

      pm=939.0d0
      pimass=140.d0
      dm=1219.d0
!      ad1=700.d0
!      ad1=700.d0	! 2H
!      ad1=774.d0	! 12C
      ad1=750.d0	! 12C
!      ad0=700.d0
      ad0=774.d0
      pi=acos(-1.d0)
      alph=1.d0/137.03604d0
      hbarc=197.32858d0
!      gamdp=110.d0
      gamdp=110.d0
!      gamsprd=140.d0
      gamsprd=20.d0
!      gamr=120.d0
      gamr=100.d0
!      gampi=5.d0
      gampi=50.d0
      qfdp=1.02d-7
      pfr=230.d0
      qmsr=4.d0*730.d0*(730.d0-390.d0)*sin(37.1d0*pi/180.d0/2.d0)**2
      qvsr=qmsr+390.d0**2
      qmsrq=4.d0*730.d0*(730.d0-115.d0)*sin(37.1d0*pi/180.d0/2.d0)**2
      qvsrq=qmsrq+115.d0**2
      na=int(a)
      if(na.eq.1)then
      qfd=qfdp
      gsprda=0.d0
      ad=ad0
!      elseif(na.lt.4)then
      elseif(na.eq.2)then
      qfd=qfdp
      gsprda=(a-1.d0)*gamsprd/3.d0
!      gsprda=gamspr
!      ad=ad0+(a-1.d0)*(ad1-ad0)/3.d0
      ad=700.d0
      else
      ad=ad1
      gsprda=gamsprd
!      qfd=qfdp
      qfd=qfdp*1.1d0	! 5/03
      endif
      thr=th*pi/180.d0
      qms=4.d0*e*(e-w)*sin(thr/2.d0)**2
      qvs=qms+w**2
      ekappa=w-qms/2.d0/pm
      cmtot2=pm**2+2.d0*pm*ekappa
c  begin delta calculation
      if(na.gt.1)then
      gamq=gamr*pf*sqrt(qvs)/pfr/sqrt(qvsrq)
      else
      gamq=0.d0
      endif
      epd=e-(dm-pm)*(dm+pm)/2.d0/pm
      epd=epd/(1.d0+2.d0*e*sin(thr/2.d0)**2/pm)
      epd=epd-epsd
      wd=e-epd
      qmspk=4.d0*e*epd*sin(thr/2.d0)**2
c     qvspk=qmspk+wd**2
c
c note width includes e-dependence,fermi broadening,& spreading
c
      wthresh=4.d0*e**2*sin(thr/2.d0)**2+pimass**2+2.d0*pimass*pm
      wthresh=wthresh/2.d0/pm
       thrshfree=1.d0+2.d0*e*sin(thr/2.d0)**2/pm
      threshd=1.d0+pf/pm+pf**2/2.d0/pm**2+2.d0*e*sin(thr/2.d0)**2/pm
      wthrfree=wthresh/thrshfree
      wthresh=wthresh/threshd
      if(w.gt.wthresh)then
        if((z.ne.1).and.(a.ne.1)) gampi = wthrfree-wthresh
      thresh=1.d0-exp( (-(w-wthresh)) / gampi )
      else
      thresh=0.d0
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
      sigd=sigd*(qms/2.d0/qvs+tan(thr/2.d0)**2)
      sigd=sigd/(qmsr/2.d0/qvsr+tan(37.1d0*pi/180.d0/2.d0)**2)
      sigd=sigd*sigmot(e,thr)/sigmot(730.d0,37.1d0*pi/180.d0)
      sigd=sigd*a                        
      sigd=sigd*thresh
      sigdel=sigd
      return
      end

*sigx
      real*8 function sigx(e,th,w,a)
      implicit none
      real*8 e,th,w,a
      real*8 alph, pi, sig0, sig1, pimass, pm, gam0, thr, qms, arg0,
     &       arg1, arg, shape, ekappa, siggam, qs, eps, flux, sigee, 
     &       fphenom, r, factor1  

      alph=1.d0/137.03604d0
      pi=acos(-1.d0)
c     sig0=111.d0*1.d-30
c modify to get agreement with Run 43204 carbon data
      sig0=100.d-4*.92
!      sig0=80.d-4
c     sig1=60.d0*1.d-27
      sig1=54.d0*1.d-1
      pimass=140.d0
      pm=939.0d0
c     gam0=550.d0
!     gam0=600.d0
      gam0=610.d0
!      gam0=700.d0
c     r=0.10d0
c     aq=250.d0
      thr=th*pi/180.d0
      if(w.lt.1.d-5)go to 4
      qms=4.d0*e*(e-w)*sin(thr/2.d0)**2
      arg0=w-qms/2.d0/pm-pimass-pimass**2/2.d0/pm
      arg1=arg0/gam0
      arg=arg1**2/2.d0
      if(arg1.gt.8.d0)then
      shape=1.d0+sig1/sig0/arg0
      elseif(arg1.lt.1.d-5)then
      shape=0.d0
      elseif(arg1.lt.0.1d0)then
      shape=sig1*arg0/2.d0/gam0**2/sig0
      else
      shape=(1.d0-exp(-arg))*(1.d0+sig1/sig0/arg0)
      endif
      ekappa=w-qms/2.d0/pm
      siggam=sig0*shape
      qs=qms+w**2
      eps=1.d0/(1.d0+2.d0*qs*tan(thr/2.d0)**2/qms)
      flux=alph*ekappa*(e-w)/2.d0/pi**2/qms/e/(1.d0-eps)
      if(flux.lt.1.d-20)flux=0.d0
      sigee=flux*siggam*fphenom(qms)**2
c     sigee=flux*siggam
!      r=0.56d0*1.6d0/(qms+pm**2)
	r=min(0.56d0*1.d6/(qms+pm**2),0.2d0)
      factor1=1.d0+eps*r
      sigee=sigee*factor1
 4    sigx=a*sigee
      return
      end

*sigr1
      real*8 function sigr1(e,th,w,a,pf,z)
      implicit none
      real*8 e,th,w,a,pf,z
      real*8 pi, pm, pimass, thr, pfr, rm, epsr, ar0, ar1, gamqfr, ar
      real*8 gamsprd, gamr, gampi, qfrp, qmsqfr, qvsqfr, qmsrr, qvsrr
      real*8 sigref, fd, sigmot, qfr, gsprda, qms, qvs, gamq, cmtot2
      real*8 wthresh, thrshfree, threshd, wthrfree, thresh, epr
      real*8 gam, sigr      
      integer na     

      pi=acos(-1.d0)
      pm=939.0d0
      pimass=140.d0
      thr=th*pi/180.d0
      pfr=230.d0
!     rm=1500.d0
      rm=1500.d0
      epsr=0.d0
      ar0=1000.d0
!      ar1=1000.d0
!      ar1=940.d0	! 2H
      ar1=1000.d0	! 12C
!      gamqfr=120.d0
      gamqfr=100.d0
!      gamsprd=140.d0
      gamsprd=0.d0
      gamr=110.d0
!      gamr=130.d0
!      gampi=5.d0
      gampi=25.d0
      qfrp=1.20d-7
      qmsqfr=4.d0*730.d0*(730.d0-115.d0)*sin(37.1d0*pi/180.d0/2.d0)**2
      qvsqfr=qmsqfr+115.d0**2
      qmsrr=4.d0*10000.d0*(10000.d0-1240.d0)*sin(6.d0*pi/180.d0/2.d0)**2
      qvsrr=qmsrr+1240.d0**2
      sigref=fd(qmsrr,ar0)**2*qvsrr
      sigref=sigref*(qmsrr/2.d0/qvsrr+tan(6.d0*pi/180.d0/2.d0)**2)
      sigref=sigref*sigmot(10000.d0,6.d0*pi/180.d0)
      na=int(a)
      if(na.eq.1)then
      qfr=qfrp
      gsprda=0.d0
      ar=ar0
!      elseif(na.lt.4)then
      elseif(na.eq.2)then
      qfr=qfrp
      gsprda=(a-1.d0)*gamsprd/3.d0
!      ar=ar0+(a-1.d0)*(ar1-ar0)/3.d0
      ar=920.d0
      else
      ar=ar1
      gsprda=gamsprd
      qfr=qfrp
      endif
      qms=4.d0*e*(e-w)*sin(thr/2.d0)**2
      qvs=qms+w**2
      if(na.gt.1)then
      gamq=gamqfr*pf*sqrt(qvs)/pfr/sqrt(qvsqfr)
      else
      gamq=0.d0
      endif
      cmtot2=pm**2+2.d0*pm*w-qms
      wthresh=4.d0*e**2*sin(thr/2.d0)**2+pimass**2+2.d0*pimass*pm
      wthresh=wthresh/2.d0/pm
       thrshfree=1.d0+2.d0*e*sin(thr/2.d0)**2/pm
      threshd=1.d0+pf/pm+pf**2/2.d0/pm**2+2.d0*e*sin(thr/2.d0)**2/pm
       wthrfree=wthresh/thrshfree
      wthresh=wthresh/threshd
      if(w.gt.wthresh)then
        if((z.ne.1).and.(a.ne.1)) gampi = wthrfree-wthresh
      thresh=1.d0-exp( (-(w-wthresh)) / gampi )
      else
      thresh=0.d0
      endif
      epr=e-(rm-pm)*(rm+pm)/2.d0/pm
      epr=epr/(1.d0+2.d0*e*sin(thr/2.d0)**2/pm)
      epr=epr-epsr
c     wr=e-epr
      gam=sqrt(gamr**2+gamq**2+gsprda**2)
      sigr=qfr*(gamr/gam)/sigref
      sigr=sigr*cmtot2*gam**2
      sigr=sigr/((cmtot2-(rm+epsr)**2)**2+cmtot2*gam**2)
      sigr=sigr*qvs*fd(qms,ar)**2
      sigr=sigr*(qms/2.d0/qvs+tan(thr/2.d0)**2)
      sigr=sigr*sigmot(e,thr)
      sigr1=a*thresh*sigr
      return
      end

*sigr2
      real*8 function sigr2(e,th,w,a,pf,z)
      implicit none
      real*8 e,th,w,a,pf,z
      real*8 pi, pm, pimass, thr, pfr, rm, epsr, ar0, ar1, gamqfr, ar
      real*8 gamsprd, gamr, gampi, qfrp, qmsqfr, qvsqfr, qmsrr, qvsrr
      real*8 sigref, fd, sigmot, qfr, gsprda, qms, qvs, gamq, cmtot2
      real*8 wthresh, thrshfree, threshd, wthrfree, thresh, epr
      real*8 gam, sigr      
      integer na

      pi=acos(-1.d0)
      pm=939.0d0
      pimass=140.d0
      thr=th*pi/180.d0
      pfr=230.d0
      rm=1700.d0
      epsr=0.d0
      ar0=1200.d0
!      ar1=1200.d0
!      ar1=1000.d0	! 2H
      ar1=1200.d0	! 12C
!      gamqfr=120.d0
      gamqfr=100.d0
!      gamsprd=140.d0
      gamsprd=0.d0
      gamr=110.d0
!      gampi=5.d0
      gampi=25.d0
      qfrp=0.68d-7
      qmsqfr=4.d0*730.d0*(730.d0-115.d0)*sin(37.1d0*pi/180.d0/2.d0)**2
      qvsqfr=qmsqfr+115.d0**2
      qmsrr=4.d0*10000.d0*(10000.d0-1520.d0)*sin(6.d0*pi/180.d0/2.d0)**2
      qvsrr=qmsrr+1520.d0**2
      sigref=fd(qmsrr,ar0)**2*qvsrr
      sigref=sigref*(qmsrr/2.d0/qvsrr+tan(6.d0*pi/180.d0/2.d0)**2)
      sigref=sigref*sigmot(10000.d0,6.d0*pi/180.d0)
      na=int(a)
      if(na.eq.1)then
      qfr=qfrp
      gsprda=0.d0
      ar=ar0
!      elseif(na.lt.4)then
      elseif(na.eq.2)then
      qfr=qfrp
      gsprda=(a-1.d0)*gamsprd/3.d0
!      ar=ar0+(a-1.d0)*(ar1-ar0)/3.d0
	 ar=980.d0
      else
      ar=ar1
      gsprda=gamsprd
      qfr=qfrp
      endif
      qms=4.d0*e*(e-w)*sin(thr/2.d0)**2
      qvs=qms+w**2
      if(na.gt.1)then
      gamq=gamqfr*pf*sqrt(qvs)/pfr/sqrt(qvsqfr)
      else
      gamq=0.d0
      endif
      cmtot2=pm**2+2.d0*pm*w-qms
      wthresh=4.d0*e**2*sin(thr/2.d0)**2+pimass**2+2.d0*pimass*pm
      wthresh=wthresh/2.d0/pm
       thrshfree=1.d0+2.d0*e*sin(thr/2.d0)**2/pm
      threshd=1.d0+pf/pm+pf**2/2.d0/pm**2+2.d0*e*sin(thr/2.d0)**2/pm
       wthrfree=wthresh/thrshfree
      wthresh=wthresh/threshd
      if(w.gt.wthresh)then
        if((z.ne.1).and.(a.ne.1)) gampi = wthrfree-wthresh
      thresh=1.d0-exp( (-(w-wthresh)) / gampi )
      else
      thresh=0.d0
      endif
      epr=e-(rm-pm)*(rm+pm)/2.d0/pm
      epr=epr/(1.d0+2.d0*e*sin(thr/2.d0)**2/pm)
      epr=epr-epsr
c     wr=e-epr
      gam=sqrt(gamr**2+gamq**2+gsprda**2)
      sigr=qfr*(gamr/gam)/sigref
      sigr=sigr*cmtot2*gam**2
      sigr=sigr/((cmtot2-(rm+epsr)**2)**2+cmtot2*gam**2)
      sigr=sigr*qvs*fd(qms,ar)**2
      sigr=sigr*(qms/2.d0/qvs+tan(thr/2.d0)**2)
      sigr=sigr*sigmot(e,thr)
      sigr2=a*thresh*sigr
      return
      end

*sig2n
      real*8 function sig2n(e,th,w,z,a,pf)
      implicit none
      real*8 e,th,w,z,a,pf
      real*8 pi, thr, dm, pimass, pm, a2, pfr, gam2n, gamqfr, gamref
      real*8 gamr, qmsr, effmass, qvsr, fd, qms, qvs, gamqf
      real*8 ekappa, cmtot2, gam, wthresh, thresh
      real*8 sig, sigref, sigkin, sigmot, sigcon

      pi=acos(-1.d0)
      thr=th*pi/180.d0
      dm=1219.d0
      pimass=140.d0
      pm=939.d0
!      a2=550.d0
      a2=500.d0
      pfr=60.d0
      gam2n=20.d0
!      gam2n=40.d0
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
      thresh=1.d0-exp( (-(w-wthresh)) / gam2n )
      else
      thresh=0.d0
      endif
      sig2n=sig2n*thresh
      return
      end

*radiate
      subroutine radiate(e,th,w,z,a,signr,sigrad,delnu)  ! 'delnu' added mar05
      implicit none
      real*8  e,th,w,z,a,signr,sigrad, delnu  ! 'delnu' added mar05
      real*8 eps,epsd,pf,spence, alph, emass, del, pi, prec, thr, arg
      real*8 qms, d1, d2, ebar, x1, x2, ans, err
      integer nsp, iflag, kf

      common /par/ eps,epsd,pf,spence
      alph=1.d0/137.03604d0
      emass=0.511d0
      del=10.d0
      pi=acos(-1.d0)
      prec=.0005d0
      thr=th*pi/180.d0
      arg=cos(thr/2.d0)**2
      spence=pi**2/6.d0-log(arg)*log(1.d0-arg)
      do 10 nsp=1,50
 10   spence=spence-arg**nsp/dble(nsp)**2   !! float changed to dble mar05
c     print 15,spence
 15   format(' spence function = ',1f14.9)
      qms=4.d0*e*(e-w)*sin(thr/2.d0)**2
      d1=(2.d0*alph/pi)*(log(qms/emass**2)-1.d0)
      d2=13.d0*(log(qms/emass**2)-1.d0)/12.d0-17.d0/36.d0
     &   -0.5d0*(pi**2/6.d0-spence)
      d2=d2*(2.d0*alph/pi)
      ebar=sqrt(e*(e-w))
      sigrad=signr*(1.d0+d2)*exp((-d1)*log(ebar/del))
c     print 20,d1,d2
 20   format(' delta1 and delta2 = ',2e14.6)
      x1=0.d0
      x2=w-del
      call rom(e,th,w,z,a,x1,x2,prec,ans,kf,delnu)  ! 'delnu' added mar05
      ans=ans*1.d-26
      iflag=0
      err=0.d0
c     print 25,x1,x2,prec,ans,err
 25   format(' x1,x2,prec,ans,err = ',5e10.3)
c     print 30,kf,iflag
 30   format(' kf,iflag = ',2i5)
      sigrad=sigrad+ans
      return
      end

*rom
      subroutine rom(e0,th0,w0,z0,a0,a,b,eps,ans,k,delnu) ! 'delnu' added mar05
      implicit none
      real*8 e0,th0,w0,z0,a0,a,b,eps,ans,delnu  ! 'delnu' added mar05
      real*8 h, fa, fb, e, x, sig, f
      integer k,m,l,iu,iv, j, j1

c  romberg method of integration
c     dimension w(50,50)
      real*8 w(50,50)
      h=b-a
      k=0
      call value(e0,th0,w0,z0,a0,a,fa,delnu)  ! 'delnu' added mar05
      call value(e0,th0,w0,z0,a0,b,fb,delnu)  ! 'delnu' added mar05
      w(1,1)=(fa+fb)*h/2.d0
    4 k=k+1
      if(k.ge.49)go to 5
      h=h/2.d0
      sig=0.d0
      m=2**(k-1)
      do 1 j=1,m
      j1=2*j-1
      x=a+dble(j1)*h                    !! float changed to dble mar05
      call value(e0,th0,w0,z0,a0,x,f,delnu)  ! 'delnu' added mar05
    1 sig=sig+f
      w(k+1,1)=w(k,1)/2.d0+h*sig
      do 2 l=1,k
      iu=k+1-l
      iv=l+1
    2 w(iu,iv)=(4.d0**(iv-1)*w(iu+1,iv-1)-w(iu,iv-1)) /
     &         (4.d0**(iv-1)-1.d0)
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
      subroutine value(e,th,w,z,a,x,f,delnu)  ! 'delnu' added mar05
      implicit none
      real*8 e,th,w,z,a,x,f, eps,epsd,pf,spence, alph, emass, pi,thr,
     &       delnu, qms1, qms2, qmsbar, f1, f2, d1, d2, ebar
      real*8 sig1,sig2, sigqfs, sigdel, sigx, sigr1, sigr2, sig2n

      common /par/ eps,epsd,pf,spence
      alph=1.d0/137.03604d0
      emass=.511d0
      pi=acos(-1.d0)
      thr=th*pi/180.d0
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
      qms1=4.d0*e*(e-x)*sin(thr/2.d0)**2
      qms2=4.d0*(e-w+x)*(e-w)*sin(thr/2.d0)**2
      qmsbar=sqrt(qms1*qms2)
      f1=(log(qms1/emass**2)-1.d0)/2.d0/(e-x)**2
      f1=f1*sig1
      f1=f1*alph*((e-x)**2+(e-w)**2)/(w-x)/pi
      f2=sig2*(log(qms2/emass**2)-1.d0)/2.d0/e**2
      f2=f2*alph*(e**2+(e-w+x)**2)/(w-x)/pi
      d1=(2.d0*alph/pi)*(log(qmsbar/emass**2)-1.d0)
      d2=13.d0*(log(qmsbar/emass**2)-1.d0)/12.d0-17.d0/36.d0
      d2=(2.d0*alph/pi)*(d2-0.5d0*(pi**2/6.d0-spence))
      ebar=sqrt((e-x)*(e-w))
      f=(f1+f2)*(1.d0+d2)*((w-x)/ebar)**d1
      return
      end
*
      subroutine get_gegm_prot(q2,ap,gepsq,gmpsq)
      implicit none
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
     1  	+2.8d0*sqrt(q2)*q2+0.84d0*q2*q2)

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
      implicit none
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
      gb = exp((-1.d0)/2.d0*((sqrt(q2)-.29d0)/.20d0)**2) 
     >     +  exp((-1.d0)/2.d0*((sqrt(q2)+.29d0)/.20d0)**2)
      ge = ge + 0.23d0*q2*gb
      gm = 1.012d0/(1.d0+q2/0.77d0)**2 - 0.012d0/(1.0d0+q2/6.8d0)**2
      gb = exp((-1.d0)/2.d0*((sqrt(q2)-.33d0)/.14d0)**2) 
     >     +  exp((-1.d0)/2.d0*((sqrt(q2)+.33d0)/.14d0)**2)
      gm = gm - 0.28d0*q2*gb
c     gd = 1.d0/(1.d0+q2/.71d0)/(1.d0+q2/.71d0)  !! not used mar,05
      q2 = 1000.d0*1000.d0*q2

c A>1 dipole correction (bound nucleons)
	ge = ge*dipcorr
	gm = gm*dipcorr

      gensq = ge**2
      gmnsq = gm**2
      endif
      return
      end
