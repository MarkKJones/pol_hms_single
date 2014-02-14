	subroutine mc_hms(p_spec, th_spec, dpp, x, y, z, dxdz, dydz, 
     >                    x_fp, dx_fp, y_fp, dy_fp, m2,ms_flag,wcs_flag, 
     >                    decay_flag, resmult, sry, srx, ok_hms)

C+_____________________________________________________________________________
!
! Monte-Carlo of HMS spectrometer.
!	Note that we only pass on real*8 variables to the subroutine.
!	This will make life easier for portability.
!
! Author: David Potterveld, March-1993
!
! Modification History:
!
!  11-Aug-1993	(D. Potterveld) Modified to use new transformation scheme in
!		which each transformation begins at the pivot.
!
!  19-AUG-1993  (D. Potterveld) Modified to use COSY INFINITY transformations.
!
!  15-SEP-1997  MODIFY stepping through spectrometer so that all drifts
!		use project.f (not transp.f), and so that project.f and
!		transp.f take care of decay.  Decay distances all assume
!		the pathlength for the CENTRAL RAY.
C-_____________________________________________________________________________

	implicit 	none

	include 'apertures.inc'
	include 'struct_hms.inc'
	include 'track.inc'
	include 'constants.inc'

*G.&I.N. STUFF - for checking dipole exit apertures and vacuum pipe to HMS hut.
	real*8 x_offset_pipes/0.0d00/,y_offset_pipes/0.0d00/
	real*8 x_offset_hut/3.5d00/,y_offset_hut/0.6d00/

! Spectrometer definitions

	integer*4 hms,sos
	parameter (hms = 1)
	parameter (sos = 2)

! Collimator (octagon) dimensions and offsets.

	real*8 h_entr,v_entr,h_exit,v_exit
	real*8 x_off,y_off,z_off

! New collimator for HMS-100 tune.
	parameter (h_entr = 4.560d00)
	parameter (v_entr = 11.646d00)
	parameter (h_exit = 4.759d00)
	parameter (v_exit = 12.114d00)

	parameter (x_off=-0.043d00)	!HMS-100 (zeroed before Fpi)
	parameter (y_off=+0.030d00)	!HMS-100 (number from gaskell)
	parameter (z_off=+40.17d00)	!HMS100 tune (dg 5/27/98)

! Math constants

	real*8 d_r,r_d,root
	parameter (d_r = pi/180.d00)
	parameter (r_d = 180.d00/pi)
	parameter (root = 0.707106781d00)		!square root of 1/2

! The arguments

	real*8 x,y,z,x0,y0			!(cm)
	real*8 dpp				!delta p/p (%)
	real*8 dxdz,dydz         		!X,Y slope in spectrometer
	real*8 x_fp,y_fp,dx_fp,dy_fp		!Focal plane values to return
	real*8 p_spec,th_spec    		!spectrometer setting
	real*8 sry				!vertical position@tgt (+y=up)
	real*8 srx				!horiz position@tgt (+x=right)
	logical ms_flag				!mult. scattering flag
	logical wcs_flag			!wire chamber smearing flag
	logical decay_flag			!check for particle decay
	logical ok_hms				!true if particle makes it

! Local declarations.

	integer*4	chan	/1/,n_classes
	logical	first_time_hms/.true./

	real*8 dpp_recon,dth_recon,dph_recon	!reconstructed quantities
	real*8 x_recon,y_recon,z_recon
	real*8 p,m2				!More kinematic variables.
	real*8 xt,yt				!temporaries
	real*8 resmult				!DC resolution factor
	real*8 xsave,ysave
	real*8 ctheta,stheta,ep
	logical dflag			!has particle decayed?
	logical ok

! Gaby's dipole shape stuff
	logical check_dipole
	external check_dipole

	integer*4 stop_id
        common /stop/ stop_id

        real*8 th_b
        logical fieldon, recon
        common /mag/ th_b, fieldon, recon

        real*8 dxdzc,dydzc
        common /fieldcheck/ dxdzc,dydzc

	real*8 bdl		!	OR - 4/04
	common/bdlsav/ bdl	!	OR - 4/04

	save		!Remember it all!

C =============================== Executable Code ============================



c	Initialize xtgt = x_tar - OR 4/04
c

c
! Initialize ok_hms to .false., reset decay flag
	ok_hms = .false.
	dflag = .false.   !particle has not decayed yet
	hSTOP_trials = hSTOP_trials + 1

! Save spectrometer coordinates and particle momentum.
        x0=x
        y0=y
        xsave=x
        ysave=y

	xs    = x
	ys    = y
	zs    = z
	dxdzs = dxdz
	dydzs = dydz
	dpps  = dpp
	p     = p_spec*(1.+dpps/100.)
c	write(*,*) '***',hSTOP_trials,x,y,z
! Read in transport coefficients.
	if (first_time_hms) then
	   call transp_init(hms,n_classes)
	   close (unit=chan)
	   if (n_classes.ne.12) stop 'MC_HMS, wrong number of transport classes'
	   first_time_hms = .false.
	endif

! Begin transporting particle.
! Do transformations, checking against apertures.


c
! Check front of fixed slit, at about 1.262+z_off meter
 
          dxdzc = dxdz
          dydzc = dydz
	    call project(xs,ys,126.2d0+z_off,decay_flag,dflag,m2,p) 
                                                            !project and decay
c
  	  if (abs(ys-y_off).gt.h_entr) then
	    hSTOP_slit_hor = hSTOP_slit_hor + 1
	    stop_id = -3
	    goto 500
	  endif
	  if (abs(xs-x_off).gt.v_entr) then
	    hSTOP_slit_vert = hSTOP_slit_vert + 1
	    stop_id = -2
	    goto 500
	  endif
	  if (abs(xs-x_off).gt. (-v_entr/h_entr*abs(ys-y_off)+3*v_entr/2)) then
	    hSTOP_slit_oct = hSTOP_slit_oct + 1
	    stop_id = -1
	    goto 500
	  endif
!Check back of fixed slit, at about 1.262+z_off+0.063 meter

	  call project(xs,ys,6.3d0,decay_flag,dflag,m2,p) !project and decay
	  if (abs(ys-y_off).gt.h_exit) then
	    hSTOP_slit_hor = hSTOP_slit_hor + 1
	    stop_id = 1
	    goto 500
	  endif
	  if (abs(xs-x_off).gt.v_exit) then
	    hSTOP_slit_vert = hSTOP_slit_vert + 1
	    stop_id = 2
	    goto 500
	  endif
	  if (abs(xs-x_off).gt. (-v_exit/h_exit*abs(ys-y_off)+3*v_exit/2)) then
	    hSTOP_slit_oct = hSTOP_slit_oct + 1
	    stop_id = 3
	    goto 500
	  endif

! Go to Q1 IN  mag bound.  Drift rather than using COSY matrices

	  call project(xs,ys,(216.075d0-126.2d0-z_off-6.3d0),decay_flag,dflag,m2,p) !project and decay
	  if ((xs*xs + ys*ys).gt.r_Q1*r_Q1) then
	    hSTOP_Q1_in = hSTOP_Q1_in + 1
	    stop_id =4
	    goto 500
	  endif

! Check aperture at 2/3 of Q1.

	  call transp(hms,2,decay_flag,dflag,m2,p,126.0d0)
	  if ((xs*xs + ys*ys).gt.r_Q1*r_Q1) then
	    hSTOP_Q1_mid = hSTOP_Q1_mid + 1
	    stop_id = 5
	    goto 500
	  endif

! Go to Q1 OUT mag boundary.

	  call transp(hms,3,decay_flag,dflag,m2,p,63.0d0)
	  if ((xs*xs + ys*ys).gt.r_Q1*r_Q1) then
	    hSTOP_Q1_out = hSTOP_Q1_out + 1
	    stop_id = 6
	    goto 500
	  endif

! Go to Q2 IN  mag bound.  Drift rather than using COSY matrices
!!	  call transp(hms,4,decay_flag,dflag,m2,p)

	  call project(xs,ys,123.15d0,decay_flag,dflag,m2,p) !project and decay
	  if ((xs*xs + ys*ys).gt.r_Q2*r_Q2) then
	    hSTOP_Q2_in = hSTOP_Q2_in + 1
	    stop_id = 7
	    goto 500
	  endif

! Check aperture at 2/3 of Q2.

	  call transp(hms,5,decay_flag,dflag,m2,p,143.67d0)
	  if ((xs*xs + ys*ys).gt.r_Q2*r_Q2) then
	    hSTOP_Q2_mid = hSTOP_Q2_mid + 1
	    stop_id = 8
	    goto 500
	  endif

! Go to Q2 OUT mag boundary.

	  call transp(hms,6,decay_flag,dflag,m2,p,71.833d0)
	  if ((xs*xs + ys*ys).gt.r_Q2*r_Q2) then
	    hSTOP_Q2_out = hSTOP_Q2_out + 1
	    stop_id = 9
	    goto 500
	  endif

! Go to Q3 IN  mag bound.  Drift rather than using COSY matrices
!!	  call transp(hms,7,decay_flag,dflag,m2,p)

	  call project(xs,ys,94.225d0,decay_flag,dflag,m2,p) !project and decay
	  if ((xs*xs + ys*ys).gt.r_Q3*r_Q3) then
	    hSTOP_Q3_in = hSTOP_Q3_in + 1
	    stop_id = 10
	    goto 500
	  endif

! Check aperture at 2/3 of Q3.

	  call transp(hms,8,decay_flag,dflag,m2,p,145.7d0)
	  if ((xs*xs + ys*ys).gt.r_Q3*r_Q3) then
	    hSTOP_Q3_mid = hSTOP_Q3_mid + 1
	    stop_id = 11
	    goto 500
	  endif

! Go to Q3 OUT mag boundary.

	  call transp(hms,9,decay_flag,dflag,m2,p,72.9d0)
	  if ((xs*xs + ys*ys).gt.r_Q3*r_Q3) then
	    hSTOP_Q3_out = hSTOP_Q3_out + 1
	    stop_id = 12
	    goto 500
	  endif

! Go to D1 IN magnetic boundary, Find intersection with rotated aperture plane.
! Aperture has elliptical form.
!!	  call transp(hms,10,decay_flag,dflag,m2,p)

	  call project(xs,ys,102.15d0,decay_flag,dflag,m2,p) !project and decay
	  call rotate_haxis(-6.0d0,xt,yt)
	  if (check_dipole(xt,yt)) then
	    hSTOP_D1_in = hSTOP_D1_in + 1
	    stop_id = 13
	    goto 500
	  endif

! Go to D1 OUT magnetic boundary.
! Find intersection with rotated aperture plane.

	  call transp(hms,11,decay_flag,dflag,m2,p,526.1d0)
c	  call rotate_haxis(6.0d0,xt,yt)
	  if (check_dipole(xs,ys)) then
	    hSTOP_D1_out = hSTOP_D1_out + 1
	    stop_id = 14
	    goto 500
	  endif

! Check a number of apertures in the vacuum pipes following the
! dipole.  First the odd piece interfacing with the dipole itself

	  if ( (((xs-x_offset_pipes)**2+(ys-y_offset_pipes)**2).gt.30.48**2)
     >           .or. (abs((ys-y_offset_pipes)).gt.20.5232) ) then
	    hSTOP_D1_out = hSTOP_D1_out + 1
	    stop_id = 15
	    goto 500
	  endif

! Check the exit of the 26.65 inch pipe

	  call project(xs,ys,64.77d0,decay_flag,dflag,m2,p) !project and decay
	  if (((xs-x_offset_pipes)**2+(ys-y_offset_pipes)**2).gt.1145.518)then
	    hSTOP_D1_out = hSTOP_D1_out + 1
	    stop_id = 16
	    goto 500
	  endif

! check exit of long (117 inch) pipe (entrance is bigger than previous pipe)
! note: Rolf claims its 117.5 but the dravings say more like 116.x
! .. so i put 117 even.  Should be a 30.62 diameter pipe

	  call project(xs,ys,298.58d0,decay_flag,dflag,m2,p) !project and decay
	  if (((xs-x_offset_hut)**2+(ys-y_offset_hut)**2).gt.1475.90)then
	    hSTOP_D1_out = hSTOP_D1_out + 1
	    stop_id = 17
	    goto 500
	  endif

! lastly check the exit of the last piece of pipe. 45.5 inches long, 30.62 dia.

	  call project(xs,ys,+114.17d0,decay_flag,dflag,m2,p) !project and decay
	  if (((xs-x_offset_hut)**2+(ys-y_offset_hut)**2).gt.2119.45)then
	    hSTOP_D1_out = hSTOP_D1_out + 1
	    stop_id = 18
	    goto 500
	  endif

! Note that we do NOT transport (project) to focal plane.  We will do this
! in mc_hms_hut.f so that it can take care of all of the decay, mult. scatt,
! and apertures.  Pass the current z position so that mc_hms_hut knows
! where to start.  Initial zposition for mc_hms_hut is -147.48 cm so that the
! sum of four drift lengths between pipe and focal plane is 625.0 cm
! (64.77+297.18+115.57+147.48=625)

! If we get this far, the particle is in the hut.

c	  hSTOP_hut = hSTOP_hut + 1
c	  stop_id = 19

! and track through the detector hut

	  call mc_hms_hut(m2,p,x_fp,dx_fp,y_fp,dy_fp,ms_flag,wcs_flag,
     >		          decay_flag,dflag,resmult,ok,-147.48d0)
	  if (.not.ok) then
	     hSTOP_hut = hSTOP_hut + 1
	     goto 500
	     stop_id = 19
	  endif
! replace xs,ys,... with 'tracked' quantities.

	  xs=x_fp
	  ys=y_fp
	  dxdzs=dx_fp
	  dydzs=dy_fp

! Reconstruct target quantities.

c         if(recon)
c     >      call mc_hms_recon(dpp_recon,dth_recon,dph_recon,x_recon,
c     >              y_recon,z_recon,sry,srx,x0,y0,p_spec,th_spec*r_d)
        if ( recon) then
c
	   call  simc_hms_recon (dpp_recon,dth_recon,dph_recon,y_recon
     >                             ,-sry,xs,dxdzs,ys,dydzs)
	   if ( fieldon) then
	      ctheta=cos(th_spec)
	      stheta=sin(th_spec)
	      ep=p_spec*(1+dpp_recon/100.)
	      call track_to_tgt(dpp_recon,y_recon,dph_recon,dth_recon,
     >          srx,sry,-ep,me,ctheta,stheta,-1,ok_hms,xs,dxdzs,ys,dydzs
     >          ,x_recon,bdl)
	      if ( .not. ok_hms) then
		 nfail_track_to_tgt =nfail_track_to_tgt + 1 
		 goto 500
              endif
	   endif
	endif
! Fill output to return to main code

	  stop_id = 20
	  dpp  = dpp_recon
	  dxdz = dph_recon
	  dydz = dth_recon
	  x    = x_recon
	  y    = y_recon
	  z    = z_recon

	  ok_hms = .true.
	  hSTOP_successes = hSTOP_successes + 1

! We are done with this event, whether GOOD or BAD.

500	continue
c
c	if (stop_id .ne. 20) write(66,*) xs,xsave,ys,ysave,z,dxdz,dydz,stop_id
c

	return
	end




