c
c	Adapted from I. Niculescu-C. Keppel-L. Stuart model
c	for inelastic deuteron cross sections
c	by O. R. 2/04
c	Fit parameters: values in I.N.'s file par_d2_jlab_slac.txt
c
!	subroutine hallc2h(e,th,w,sigdel,sigr1,sigr2,sigx)
	subroutine hallc2h(e,th,enu,sigdel,sigr1,sigr2,sigx)
c	implicit real*8 (a-h,o-z)
	implicit none
        real*8 e,th,enu,sigdel,sigr1,sigr2,sigx
	real*8 wr(3), gr(3), cr(3,4), cnr(3,5), qk(3)
	real*8 bw(3), polyr(3), polynr(3)
        real*8 pi, pm, alpha, wpi, thr, sinth22, pm2, ef,
     &         qsq, ww, w, k, eps, g, sigmanrt, sigmanrl, dipole2
        integer j, jp

	parameter (pi=3.1415928d0)
	parameter (pm = 939.d0)	
	parameter (alpha = 1.d0/137.036d0)
!	parameter (wpi = 1.078d0)
c       parameter (wpi = 0.982d0)       ! wpi with kF = 46 MeV/c              
	parameter (wpi = 1.038d0)  !! <== check modified wpi, 2006 nov

	data wr/1.210098d0, 1.509d0, 1.675d0/

	data gr/0.16559d0, 0.133d0, 0.191d0/

	data qk/0.0d0, -1.33d-3, 9.04d-3/

	data cr/ 88.31d0,   17.276d0,  10.571d0,
     1  	 179.83470d0,  14.539d0,  17.876d0,
     1  	 -49.3425d0,   0.0d0,     0.0d0,                 
     1  	   4.419d0,   0.0d0,     0.0d0/
	
	data cnr/-1596.6d0,  3161.6d0, -1294.7d0,
     1	          3765.5d0, -5091.2d0,  1484.7d0,
     1  	  1854.5d0,  1327.6d0,  -905.3d0,
     1	          -203.6d0,  -248.9d0,   285.3d0,
     1  	     1.67d0,   26.0d0,   -24.5d0/


	thr = th*pi/180.d0                  
	sinth22 = (sin(thr/2.d0))**2
	pm2 = pm*pm
	ef = e - enu
	qsq = 4.d0*e*ef*(sin(thr/2.d0))**2
	ww = pm**2 + 2.d0*pm*enu - qsq
!	ww = w*w
!	ef = (pm2 + 2.d0*pm*e - ww)/(2.d0*pm + 4.d0*e*sinth22)
!	qsq = 4.d0*e*ef*sinth22
!	enu = e - ef
	k = (ww-pm2)/(2.d0*pm)

	if ( ww .le. 0) then
 	 sigdel = 0.d0
	 sigr1 = 0.d0
	 sigr2 = 0.d0
	 sigx = 0.d0
	 return
        else
	   w = sqrt(ww)*1.d-3	! MeV -> GeV
 	   if(w.le.wpi) then
 	    sigdel = 0.d0
	    sigr1 = 0.d0
	    sigr2 = 0.d0
	    sigx = 0.d0
	    return
	    end if
	endif

	eps = 1.d0/(1.d0+2.d0*(enu**2/qsq+1.d0)*(tan(thr/2.d0))**2)
	g = alpha*k/(2.d0*pi**2*qsq)*(ef/e)/(1.d0 - eps)
	qsq =qsq*1.d-6	! MeV**2 -> GeV^2


	sigmanrt = 0.d0
	sigmanrl = 0.d0


	do j = 1, 3

	bw(j) = gr(j)/((w-wr(j)*(1.d0+qk(j)*qsq))**2+gr(j)**2/4.d0)

	polyr(j) = 0.d0

	do jp = 1, 4
	polyr(j) = cr(j,jp)*qsq**(jp-1)+polyr(j)
	end do

	polynr(j) = 0.d0

	do jp = 1, 5
	polynr(j) = cnr(j,jp)*qsq**(jp-1)+polynr(j)
	end do

	sigmanrt = polynr(j)*sqrt((w-wpi)**(2*j-1)) + sigmanrt

	end do

	dipole2 = 1.0d0/(1.d0+qsq/0.71d0)**4

	sigmanrl = sigmanrt*0.25d0/sqrt(qsq)

	sigdel = g*bw(1)*polyr(1)*dipole2*1.d-4	!1.d-4 to combine with scale in main for nb/sr/MeV
	sigr1 = g*bw(2)*polyr(2)*dipole2*1.d-4
	sigr2 = g*bw(3)*polyr(3)*dipole2*1.d-4
	sigx = g*(sigmanrt+eps*sigmanrl)*dipole2*1.d-4
	return
	end
