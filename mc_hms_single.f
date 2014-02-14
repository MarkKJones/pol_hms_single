C	program mc_hms_single
C------------------------------------------------------------------------------
c Monte-Carlo of HMS spectrometer using uniform illumination.
c   This version uses TRANSPORT right-handed coordinate system.
c
c Author: David Potterveld, March-1993
c
c Modification History:
c
c  11-Aug-1993	(D. Potterveld) Modified to use new transformation scheme in
c		which each transformation begins at the pivot.
c
c  19-AUG-1993  (D. Potterveld) Modified to use COSY INFINITY transformations.
C------------------------------------------------------------------------------

	implicit none

	include 'struct_hms.inc'
	include 'constants.inc'
        include 'track.inc'
	include 'radc.inc'

	integer*4 pawc_size
	parameter (pawc_size = 4000000)
	common /pawc/ hbdata(pawc_size)
	integer*4 hbdata
	common /QUEST/IQUEST(100)
	integer*4 IQUEST
	integer nvars     
!	parameter (nvars=28)
	parameter (nvars=29)
	character*8 hut_nt_names(nvars)/
     >	'numer_wt', 'A',        'W',        'hsdeltai', 'hsdelta',
     >  'hsthetai', 'hstheta',  'thi',      'th',       'hsphi_i',
     >  'hsphi',       'hsxptari', 'hsxptar',  'hsyptari', 'hsyptar',  
     >  'Q2',       'hsxfp',    'hsyfp',    'hsxpfp',   'hsypfp',   
     >  'hsxtari',  'hsxtar',   'hsytari',  'hsytar',   'hszbeam',
     >  'zbeami',   'ybeami',   'xbeami',   'elev'/
	real*4 hut(nvars) ! needs to be real*4 for CERN compatiblity

C command line arguments

	integer iargc,lnblnk	!intrinsic functions
        character*80 rawname,inpname,tarname,hbname,outname,logname,
     >    hbname0, ntp_number_char
        integer ntp_number, ntp_events, ntp_events_max
        parameter(ntp_events_max=500000)     !! good up to 500,000.

C target material

	integer*4 i,j,k,l
	integer*4 nmat !number of materials comprising target
	integer*4 imat !material number from which scattering occurs
	real*8 m_n
	real*8 zz(20),a(20),thick(20),density(20),radlen(20),tar_zoff(20) 
	real*8 current,lumin(20),prob(20),int(20)
	real*8 ver(20) !thickness*density for each target material
	parameter(current=0.100)   !100 nA

C from input file

	character*70 str_line
	logical*4 iss, rd_int, rd_real
	integer*4 last_char, tmp_int, chanin /1/, chanout /2/
        integer*4 n_trials, trial
	real*8 e,p_spec,th_spec,cos_ts,sin_ts
	real*8 ep_min, ep_max
        real*8 dpp_lo,dpp_hi !mc limits on dpp
	real*8 th_lo,th_hi,ph_lo,ph_hi !mc limits on th, ph, stored in radians
        real*8 raster_radius,raster_xoff,raster_yoff !in cm
	real*8 cut_dpp,cut_dth,cut_dph,cut_z    
	integer*4 p_flag
        logical ms_flag /.false./ 
        logical wcs_flag /.false./
	logical*4 decay_flag /.false./

C event and cross section

	real*8 r,phi,x,y,z,dpp,dxdz,dydz,dydz_off,ue(3) 	
	real*8 dpp_init
        real*8 dth_init,dph_init !init angles in hms coord
	real*8 xtar_init,ytar_init,ztar_init !scat point in hms coord
	real*8 th_ev !full lab scattering angle, degrees
	real*8 th    !full lab scattering angle, degrees
	real*8 theta_ev !in plane scattering angle
	real*8 phi_ev   !out of plane scattering angle
	real*8 hsphi_init,hsphi ! azimuthal angle
	real*8 th_ev_init,th_ev_recon,cos_ev,cos_ev_recon
        real*8 theta_ev_init,theta_ev_recon,phi_ev_init,phi_ev_recon
	real*8 before(20),after(20) !thickness, by matl, bef & after scatter
        real*8 tb,ta !total rad length before and after scattering point
	real*8 enu,ep,Q2_vertex,W_vertex,erad,eprad    
	real*8 m2,sigrad,domega, denergy, normrate, sigma_qfs
	real*8 rate(20),trials(20),sigma_hms(20),rate_hms,ave_rate
                           !rates and trials for each target layer
	real*8 delnu
        parameter (delnu=14.d0)   !FWHM of proton elastic peak in MeV 

C initial and reconstructed track quantities

	real*8 resmult
	real*8 dpp_recon,dth_recon,dph_recon,xtar_recon,ytar_recon,ztar_recon
	real*8 w_recon,q2_recon,ep_recon
	real*8 x_fp,y_fp,dx_fp,dy_fp		
	real*8 dpp_var(2),dth_var(2),dph_var(2),ztg_var(2)
        real*8 t1,t2,t3,t4
	real*8 xbeam_init,ybeam_init,zbeam_init !scatt pos in beam coord
	logical*4 ok_spec	         

	real*8 theta_pol,escat

	real*8 cer_eff,dc_eff
C miscellaneous

	integer prev_suc
	real*8 zrand
	double precision grnd

	common /stop/ stop_id
        integer*4 stop_id

        real*8 th_b
        logical fieldon /.false./
        logical recon /.false./
        common /mag/ th_b, fieldon, recon

        real*8 dxdzc,dydzc
        common /fieldcheck/ dxdzc,dydzc

	logical*4 radiate /.false./      !radiate simc way
        logical*4 rad_qfs /.false./      !QFS internal radiation
        real*8 rad_weight
        logical success

	real*8 zHMS	! for tar_zoff.ne.0 OR 9/03
!	real*8 dysdzs,dxsdzs	! angles in HMS system rotated from BEAM 9/03

         character*17 trg_field_map /'trg_field_map.dat'/
         real*8 btheta_diff

	real*8 bdl!	OR - 4/04
	common/bdlsav/ bdl!	OR - 4/04
!
	real*8 proton_wgt,elastic_event
	real*8 enu_vertex,ebeam_vertex,ebeam,eprime_vertex
	real*8 etemp,eptemp
	real*8 ep_rad_max,ep_rad_min
	real*8 scalein
	real*8 packfrac
c
	real*8 hscat_win_den,hscat_win_z,hscat_win_a,hscat_win_len
	real*8 hscat_win_eloss
	real*8 temp_eloss
	real*8 tot_ran_ebeam_loss,tot_mp_ebeam_loss
	real*8 ebeam_ran,ebeam_mp
	real*8 p_temp
	real*8 targ_win_eloss,targ_he_eloss,front_eloss
	integer c_in_tar,he_in_tar,tartype
	real*8 tar_a(6),tar_z(6),tar_t(6),tar_den(6)
c	real*8 tar_lrad(6)
	real*8 tar_len(6)
	real*8 tot_ran_escat_eloss
	real*8 path,x_rear_wall,x_front_wall,z_side_wall
	real*8 path_in_he,path_in_he_after,path_tail,path_4k
	real*8 mscat_len
	real*8 beta_temp,gamma_temp,velocity,th_temp
	real*8 tot_mp_escat_eloss
	real*8 rear_wall_zpos,front_wall_zpos
	data tar_a / 12.0d0,12.0d0,4.0d0,27.0d0,18.04d0,21.04d0/
	data tar_z /6.0d0,6.0d0,2.0d0,13.0d0,10.0d0,10.0d0/
c	data tar_lrad /3.56d0,3.56d0,0.615d0,0.2d0,3.18d0,3.14d0/ ! 
	data tar_t /1.518d0,1.518d0,0.58d0,.048d0,1.377d0,1.585d0/ ! g/cm^2
	data tar_den /2.2d0,2.2d0,0.145d0,2.7d0,.917d0,1.056d0/  ! g/cm^3
c
	real*8 sig_pb,normrate2
	integer itime
        character	timestring*30
        double precision genrand_real1
c
	save		!Remember it all
c
          itime=time8()
c   	  call ctime(itime,timestring)
	  call sgrnd(itime)

C----------------------------- Executable Code -------------------------------

c
	do i=1,6
	   tar_len(i) = tar_t(i)/tar_den(i)
	enddo
C initialize, where do we loose events?
	

	do i=1,21  !number of elem in hstop array
	  hSTOP(i)=0
	end do

	stop_id = 0 

	do i=1, 20
          rate(i) = 0
	end do

C get the command line arguments

	if(iargc().ne.3) goto 9999

	call getarg(1,rawname)
	i=lnblnk(rawname)
        inpname = 'infiles/'//rawname(1:i)//'.inp'
	outname = rawname(1:i)//'_'
	
	call getarg(2,rawname)
	i=lnblnk(rawname)
        tarname = 'tarfiles/'//rawname(1:i)//'.tgt'
	j=lnblnk(outname)
	outname = outname(1:j)//rawname(1:i)//'_'

	call getarg(3,rawname)
	i=lnblnk(rawname)
	j=lnblnk(outname)
c
c       Let's create multiple ntuple files for each MC run.
c
	hbname0 = 'outfiles/hbook/'//outname(1:j)//rawname(1:i)
        k=lnblnk(hbname0)
        ntp_number=1               !! First ntuple file number
        call DIGIT2STRING(ntp_number, ntp_number_char)
        l=lnblnk(ntp_number_char)
	hbname = hbname0(1:k)//'.'//ntp_number_char(1:l)//'.hbook'
	logname = 'outfiles/'//outname(1:j)//rawname(1:i)//'.out'

c       print *, 'hbname0=', hbname0
c       print *, 'hbname=', hbname
c       stop



C initialize HBOOK/NTUPLE if necessary

	write(*,*) ' open hbook',hbname
	call hlimit(pawc_size)
	IQUEST(10) = 64000
	call hropen(30,'HUT',hbname,'NQ',1024,i)
	write(*,*) ' status = ',i
	if (i.ne.0) then
	  print *,'HROPEN error: istat = ',i
	  stop
        endif
	write(*,*) ' nvars',nvars
        call hbookn(1,'HUT NTUPLE',nvars,'HUT',1024,hut_nt_names)	  

C open Output file
	write(*,*) ' open output ',logname

	open (unit=chanout,status='replace',file=logname)

C read input file
	write(*,*) ' opening file',inpname
	open(unit=chanin,status='old',file=inpname)

C read in real*8's from setup file, strip off comment lines begin with '!'

	str_line = '!'
	do while (str_line(1:1).eq.'!')
	  read (chanin,1001) str_line
	enddo
	print *,'str=',str_line,'<'

! total event throw

	iss=rd_int(str_line,n_trials)
	print *,iss,n_trials,'>',str_line,'<'
	if(.not. iss) stop 'ERROR (ntrials) in setup!'
	print *,str_line(1:last_char(str_line))

! beam energy: e

	read (chanin,1001) str_line
	print *,str_line(1:last_char(str_line))
	iss = rd_real(str_line,e)
	if (.not.iss) stop 'ERROR (Beam energy) in setup!'

! spectrometer momentum

	read (chanin,1001) str_line
	print *,str_line(1:last_char(str_line))
	iss = rd_real(str_line,p_spec)
	if (.not.iss) stop 'ERROR (Spec momentum) in setup!'

! electron spectrometer angle

	read (chanin,1001) str_line
	print *,str_line(1:last_char(str_line))
	iss = rd_real(str_line,th_spec)
	if (.not.iss) stop 'ERROR (e- Spec theta) in setup!'
	th_spec = abs(th_spec)/degrad
	cos_ts = cos(th_spec)
	sin_ts = sin(th_spec)

! MC limits (half width's for dp,th,ph)

	read (chanin,1001) str_line
	iss = rd_real(str_line,dpp_lo)
	if (.not. iss) stop 'ERROR dp/p lower limit in setup!'

	read (chanin,1001) str_line
	iss = rd_real(str_line,dpp_hi)
	if (.not. iss) stop 'ERROR dp/p upper limit in setup!'
	
	read (chanin,1001) str_line
	iss = rd_real(str_line,th_lo)
	if (.not. iss) stop 'ERROR theta lower limit in setup!'
	th_lo=th_lo/1000. !convert from mr to radians

	read (chanin,1001) str_line
	iss = rd_real(str_line,th_hi)
	if (.not. iss) stop 'ERROR theta upper limit in setup!'
	th_hi=th_hi/1000. !convert from mr to radians
	
	read (chanin,1001) str_line
	iss = rd_real(str_line,ph_lo)
	if (.not. iss) stop 'ERROR phi lower limit in setup!'
	ph_lo=ph_lo/1000. !convert from mr to radians

	read (chanin,1001) str_line
	iss = rd_real(str_line,ph_hi)
	if (.not. iss) stop 'ERROR phi upper limit in setup!'
	ph_hi=ph_hi/1000. !convert from mr to radians

! Raster position and size
	read (chanin,1001) str_line
	iss=rd_real(str_line,raster_radius)
	if(.not. iss) stop 'ERROR Raster radius in setup!'

	read (chanin,1001) str_line
	iss=rd_real(str_line,raster_xoff)
	if(.not. iss) stop 'ERROR Raster x offset in setup!'

	read (chanin,1001) str_line
	iss=rd_real(str_line,raster_yoff)
	if(.not. iss) stop 'ERROR Raster y offset in setup!'

! solid angle


! cuts on reconstructed quantities

	read (chanin,1001) str_line
	print *,str_line(1:last_char(str_line))
	if (.not.rd_real(str_line,cut_dpp)) stop 'ERROR (CUT_DPP) in setup!'

	read (chanin,1001) str_line
	print *,str_line(1:last_char(str_line))
	if (.not.rd_real(str_line,cut_dth)) stop 'ERROR (CUT_DTH) in setup!'

	read (chanin,1001) str_line
	print *,str_line(1:last_char(str_line))
	if (.not.rd_real(str_line,cut_dph)) stop 'ERROR (CUT_DPH) in setup!'

	read (chanin,1001) str_line
	print *,str_line(1:last_char(str_line))
	if (.not.rd_real(str_line,cut_z)) stop 'ERROR (CUT_Z) in setup!'

! read in flag for particle type

	read (chanin,1001) str_line
	print *,str_line(1:last_char(str_line))
	if (.not.rd_int(str_line,p_flag)) stop 'ERROR: p_flag in setup!'

! read in flag for multiple scattering

	read (chanin,1001) str_line
	print *,str_line(1:last_char(str_line))
	if (.not.rd_int(str_line,tmp_int)) stop 'ERROR: ms_flag in setup!'
	if (tmp_int.eq.1) ms_flag = .true.

! read in flag for wire chamber smearing

	read (chanin,1001) str_line
	print *,str_line(1:last_char(str_line))
	if (.not.rd_int(str_line,tmp_int)) stop 'ERROR: wcs_flag in setup!'
	if (tmp_int.eq.1) wcs_flag = .true.

! read in QFS internal radiation flag

	read (chanin,1001) str_line
	print *,str_line(1:last_char(str_line))
	if (.not.rd_int(str_line,tmp_int)) stop 'ERROR: QFS internal radiation flag in setup!'
	if (tmp_int.eq.1) rad_qfs = .true.

! read in target magnetic field switch

	read (chanin,1001) str_line
	print *,str_line(1:last_char(str_line))
	if (.not.rd_int(str_line,tmp_int)) stop 'ERROR: target field switch in setup!'
	if (tmp_int.eq.1) fieldon = .true.

! read in target magnetic field angle

	read (chanin,1001) str_line
	print *,str_line(1:last_char(str_line))
	iss = rd_real(str_line,th_b)
	if (.not.iss) stop 'ERROR target field angle in setup!'
!fixme: target field angle is not converted to radians.  is this OK?

! read in reconstruction flag

	read (chanin,1001) str_line
	print *,str_line(1:last_char(str_line))
	if (.not.rd_int(str_line,tmp_int)) stop 'ERROR: recon flag in setup!'
	if (tmp_int.eq.1) recon = .true.

! read in radiation flag (simc way)

	read (chanin,1001) str_line
	print *,str_line(1:last_char(str_line))
	if (.not.rd_int(str_line,tmp_int)) stop 'ERROR: radiation flag in setup!'
	if (tmp_int.eq.1) radiate = .true.

! read in scale factor for A>2 inelastic cross section.
	read (chanin,1001) str_line
	print *,str_line(1:last_char(str_line))
	iss = rd_real(str_line,scalein)
	if (.not.iss) stop 'ERROR: A>2 inelastic scale in setup!'

! read in packing fraction ( absolute)
	read (chanin,1001) str_line
	print *,str_line(1:last_char(str_line))
	iss = rd_real(str_line,packfrac)
	if (.not.iss) stop 'ERROR:  packing fraction in setup!'

C read target material file ( assume that it is for 50% packing fraction)

	open(unit=12, file=tarname, status='old')
	read(12,*) nmat
	do i=1, nmat
	  read(12,*) zz(i), a(i), thick(i), density(i), radlen(i), tar_zoff(i)
	  if (i .eq. 1 .and. packfrac .gt. 0) density(i)=density(i)*(packfrac/0.5)
	  if (i .eq. 2 .and. packfrac .gt. 0) density(i)=density(i)*(packfrac/0.5)
	  if (i .eq. 3 .and. packfrac .gt. 0) density(i)=density(i)*((1-packfrac)/0.5)
	end do
	close(12)

	do i=1, nmat
          lumin(i) = thick(i)*density(i)*current/a(i)*N_A/Q_E*1000. !per nbarn 
	  prob(i) = density(i)*thick(i)
          write(*,*) 'lumin (per nbarn for 1uC) = ',a(i),lumin(i)*10  
	end do

	do i=1, nmat
	  int(i) = 0
	  do j=1, i
            int(i) = int(i) + prob(j)
	  end do
        end do 

	do i=1, nmat
	   int(i)=int(i)/int(nmat)   
        end do 

	ver(1) = int(1)
	do i=2, nmat 
	   ver(i)=int(i)-int(i-1)   ! ver(i)=rho(i)*t(i)/sum_i[rho(i)*t(i)]
        enddo
c
c determine target type for energy loss
c       1 = c + he
c       2 = c
c       3 = MT + He
c       4 = MT + no He
c       5 = Nh3
c       6 = Nd3
c
	c_in_tar = 0
	he_in_tar = 0
	tartype=0
	do i=1,nmat
          if ( zz(i) .eq. 1 .and. a(i) .eq. 1) tartype=5 
          if ( zz(i) .eq. 1 .and. a(i) .eq. 2) tartype=6 
          if ( zz(i) .eq. 6 .and. a(i) .eq. 12) c_in_tar=1 
          if ( zz(i) .eq. 2 .and. a(i) .eq. 4) he_in_tar=1 
	enddo
        if (he_in_tar .eq. 1 .and. c_in_tar .eq. 1 ) tartype=1
        if (he_in_tar .eq. 0 .and. c_in_tar .eq. 1 ) tartype=2
        if ( tartype .eq. 0 .and. he_in_tar .eq. 1 ) tartype=3
        if ( tartype .eq. 0 .and. he_in_tar .eq. 0 ) tartype=4
	if ( tartype .eq. 0) then
	   write(*,*) ' Cannot determine target type'
	   stop
	endif
	write(*,*) ' target type = ',tartype
c

C set particle masses

	m2 = me2			!default to electron
	if(p_flag.eq.0) then
	  m2 = me2
	else if(p_flag.eq.1) then
	  m2 = mp2
	else if(p_flag.eq.2) then
	  m2 = md2
	else if(p_flag.eq.3) then
	  m2 = mpi2
	else if(p_flag.eq.4) then
	  m2 = mk2
	else
	  stop 'ERROR: bad particle type in input file'
	endif
c
c initial target field map
	if ( fieldon) then
	write(*,*) 'Opening  target file map = ',trg_field_map
	btheta_diff = th_b - th_spec*degrad
	call trginit(trg_field_map,btheta_diff,0.d00,0.d00,0.d00)
	endif
c
c scattering chamber window
      	 hscat_win_den= 2.70    ! grams/cm**3
      	 hscat_win_z= 13.0
      	 hscat_win_a= 27.0
      	 hscat_win_len= (0.008+0.0015+0.0015+0.002)*2.54 ! RSS para
c

c
C----------------------------------------------------------------------------C
C                           Top of Monte-Carlo loop                          C
C----------------------------------------------------------------------------C
c
	nfail_track_from_tgt = 0
	ntp_events = 0
	prev_suc = 1
c
	do trial = 1,n_trials

	ep_min = p_spec*(1.+0.01*dpp_lo)
	ep_max = p_spec*(1.+0.01*dpp_hi)
	domega = (th_hi-th_lo)*(ph_hi-ph_lo)
	denergy = ep_max-ep_min

	ep_rad_min = p_spec*(1.-0.305)
	ep_rad_max = p_spec*(1.+0.305)
 	  if(mod(trial,1000).eq.0) then 
	    write(*,*)'ev= ',trial,' ,success= ',hSTOP_successes, ' ,av. rate= ', ave_rate
	  
C check if the number of successes froze. If this is the case, stop

     	    if (hSTOP_successes.ge.9000 .and. hSTOP_successes.eq.prev_suc) then 
	      print *, 'number of successes froze!'
	      goto 533
	    end if
	    prev_suc = hSTOP_successes
	  end if


!choose the x,y,z (in cm) location of the interaction point

	  !pick rnd point in our circular raster pattern for x and y pos
          r = grnd()*raster_radius
          phi = grnd()*2.*pi
	  x = raster_xoff + sqrt(r)*cos(phi)
	  y = raster_yoff + sqrt(r)*sin(phi)

	  !pick rnd target material, weighted by density*length, store in imat
	  zrand = grnd()
	  do i=1, nmat
	    if (zrand.le.int(i)) then
              imat = i
	      goto 99
	    end if
	  end do
 99	  continue
	  !pick rnd point inside this matl for z pos of interaction point
	  z = tar_zoff(imat) + (grnd()-0.5)*thick(imat)
c
c
c
c to mimic the ENGINE for correction to energy loss
          beta_temp  = 1./sqrt(1.+(.511d0/e)**2)
          gamma_temp = sqrt(1.+(.511d0/e)**2)/(.511d0/e)
          velocity   = log(beta_temp*gamma_temp)/log(10.)
	  th_temp=hscat_win_len*hscat_win_den
	  call loss(0,hscat_win_z,hscat_win_a,th_temp,hscat_win_den,velocity,hscat_win_eloss)
	  tot_mp_ebeam_loss =hscat_win_eloss 
c
c
	  if ( tartype .eq. 5 .or. tartype .eq. 6) then
	  th_temp=.00381*2.7
	  call loss(0,13.d0,26.d0,th_temp,2.7d0,velocity,targ_win_eloss)
	  tot_mp_ebeam_loss = tot_mp_ebeam_loss + targ_win_eloss
	  th_temp=1.25*.145
	  call loss(0,2.d0,4.d0,th_temp,.145d0,velocity,targ_he_eloss)
	  tot_mp_ebeam_loss = tot_mp_ebeam_loss + targ_he_eloss
	  elseif ( tartype .eq. 1) then ! C + He  len  = (4 - .8)/2
	  th_temp=1.6d0*0.145d0
	  call loss(0,2.d0,4.d0,th_temp,.145d0,velocity,targ_he_eloss)
	  tot_mp_ebeam_loss = tot_mp_ebeam_loss + targ_he_eloss
	  endif
c
	  th_temp=tar_len(tartype)*tar_den(tartype)/2.
	  call loss(0,tar_z(tartype),tar_a(tartype),th_temp,tar_den(tartype),velocity,front_eloss)
	  tot_mp_ebeam_loss = tot_mp_ebeam_loss + front_eloss
c
	  ebeam_mp = e - tot_mp_ebeam_loss

c
c determine random energy loss for beam
c
	  call enerloss_new(hscat_win_len,hscat_win_den,
     >hscat_win_z,hscat_win_a,e,.511d0,1,hscat_win_eloss)
	  tot_ran_ebeam_loss =hscat_win_eloss 
c
	  do i=1,nmat
	     temp_eloss=0
             if (tar_zoff(i)-thick(i)/2.ge.z) then
	      before(i)=0.
	     else if (tar_zoff(i)+thick(i)/2.le.z) then 
	      before(i)=thick(i)
	    else 
	      before(i)=z-(tar_zoff(i)-thick(i)/2.) 
	    end if
	    if ( before(i) .ne. 0 ) then
               call enerloss_new(before(i),density(i),
     >zz(i),a(i),e,.511d0,1,temp_eloss)
	  tot_ran_ebeam_loss = tot_ran_ebeam_loss + temp_eloss
	    endif
          enddo
	  ebeam_ran = e -  tot_ran_ebeam_loss
c

C pick scattering angles and dpp from independent, uniform distributions.
C dxdz and dydz in HMS TRANSPORT coordinates

!	  dydz_off = tar_zoff(imat)*sin(th_spec)/166.37
	  zHMS = 166.37		! correct horiz. angle offset for zoff OR 9/03
				! = 1.26m + zoff (~40cm)
	  dydz_off = tar_zoff(imat)*sin(th_spec)/(zHMS-tar_zoff(imat)*cos(th_spec)) ! 9/03
c	  if ( dydz_off .ne. 0) write(*,*) ' imat = ',imat,'dydz_off = ',dydz_off
	  dydz = dydz_off + grnd()*(th_hi-th_lo)+th_lo
	  dxdz =            grnd()*(ph_hi-ph_lo)+ph_lo
c mkj Jan/13/2005 for a=1 have half events elastic and half inelastic
	  proton_wgt = 1.
	  elastic_event=0.
          if ( a(imat) .eq. 1 .and. zz(imat) .eq. 1 ) then ! proton target
          if ( grnd() .le. 0.5 ) then 
	    dpp  = grnd()*(dpp_hi - dpp_lo) + dpp_lo
	    proton_wgt = 2.
	  else
            theta_pol = acos( (cos_ts + dydz*sin_ts)
     +                        / sqrt( 1. + dxdz**2 + dydz**2 ) )
	    escat = 938.27*ebeam_ran/(ebeam_ran*(1-cos(theta_pol))+938.27)
	    dpp = (escat-p_spec)/p_spec*100.
	    proton_wgt = 2.  
	    elastic_event=1.0
          endif
          else
	    dpp  = grnd()*(dpp_hi - dpp_lo) + dpp_lo
          endif

		
C unit vector along scattered electron (in HMS system)
c   mkj change because radiative code wants unit
c       vector relative to beam coordinate system
c          ue(1) = dxdz/sqrt(1.0+dxdz**2+dydz**2)   !DOWN
c          ue(2) = dydz/sqrt(1.0+dxdz**2+dydz**2)   !LEFT
c          ue(3) = 1.0/sqrt(1.0+dxdz**2+dydz**2)    !DOWNSTREAM
	ue(1) = dxdz/sqrt(1.0+dxdz**2+dydz**2)
	ue(2) = (dydz*cos_ts-sin_ts)/sqrt(1.0+dxdz**2+dydz**2)
	ue(3) = (dydz*sin_ts+cos_ts)/sqrt(1.0+dxdz**2+dydz**2)
	
C transform from target to HMS (TRANSPORT) coordinates: xs(down), ys(HMS left) 
C and zs(hms forward). Note that this assumes that HMS is on the right-hand 
C side of the beam line (looking downstream)
c
c  The beam coordinate system is +y (vertical up) , +z ( downstream towards dump)
c    and +x is beam right
c
	  xs = -y
	  ys = -(x * cos_ts) + z * sin_ts
	  zs = z * cos_ts + x * sin_ts
C version for spectrometer on the left-hand side:
!	  xs = -y
!	  ys = x * cos_ts - z * sin_ts
!	  zs = z * cos_ts + x * sin_ts
	  dpps  = dpp
	  dxdzs = dxdz
	  dydzs = dydz


C scattering angle
	  cos_ev   = (cos_ts+dydzs*sin_ts)/sqrt(1+dydzs**2+dxdzs**2)
	  th_ev    = acos(cos_ev)   !Lab scattering angle          
	  theta_ev = acos((cos_ts+dydzs*sin_ts)/sqrt(1+dydzs**2))  !inplane ang
          phi_ev   = asin(dxdzs/sqrt(1+dxdzs**2+dydzs**2))   !outplane angle
	  hsphi_init = atan2((dydzs*cos_ts-sin_ts),dxdzs)
           hsphi_init =  hsphi_init + 270.*3.14159/180.
 
!	print *, th_ev*degrad,theta_ev*degrad,phi_ev*degrad	!9/03
          th_ev_init    = th_ev
          theta_ev_init = theta_ev
          phi_ev_init   = phi_ev
	  ep = p_spec*(1+0.01*dpps)
C determine rad. length before and after vertex 

	  tb = 0
	  ta = 0 
	  tot_ran_escat_eloss = 0.0
	  mscat_len = 0
	  do i=1, nmat
	    if (tar_zoff(i)-thick(i)/2.ge.z) then
	      before(i)=0.
	    else if (tar_zoff(i)+thick(i)/2.le.z) then 
	      before(i)=thick(i)
	    else 
	      before(i)=z-(tar_zoff(i)-thick(i)/2.) 
	    end if
	    after(i)=thick(i)-before(i)	     
	    tb = tb + before(i)*density(i)/radlen(i)   !radlen(i) in g/cm**2
	    ta = ta + after(i)*density(i)/radlen(i)
            if (after(i) .ne. 0) then
	       after(i) = after(i)/cos_ev
             call enerloss_new(after(i),density(i),
     >zz(i),a(i),ep,.511d0,1,temp_eloss) ! He before target
              tot_ran_escat_eloss  = tot_ran_escat_eloss + temp_eloss
            endif
	  end do
c add rad len of beam exit, 5cm of air, OVC ( for perp it is Be, for para it is Al) 
c       and LN2 shield to the before rad length
          if ( th_b .eq. -90) then
	  tb = tb + (.015+.015)*2.54*1.85/64.19 + 5*.0013/36.66 + .0015*2.54*2.7/24.011 
          else 
	  tb = tb + .015*2.54*1.85/64.19 + 5*.0013/36.66 + (.008+.0015)*2.54*2.7/24.011 
	  endif
c
	  ta = ta/cos_ev	!'before' is along beam axis, 'after' at angle to it	  
	  mscat_len = ta
c add rad len of LN2 sheild, OVC, 10cm of air, HMS kevlar and mylar
          if ( th_b .eq. -90) then
	  ta = ta +  (.0015+.016)*2.54*2.7/24.011 + 10.*.0013/36.66 + .038*.74/55.2 + .0127*1.39/39.95
	  else
	  ta = ta +  (.0015+.008)*2.54*2.7/24.011 + 10.*.0013/36.66 + .038*.74/55.2 + .0127*1.39/39.95
	  endif
c
c  random energy loss for scattered electron in NH3 or ND3 is special case
               if (tartype .eq. 5 .or. tartype .eq. 6) then
         	  tot_ran_escat_eloss = 0.0
	           mscat_len = 0
		  path = 0.0
		    z_side_wall = -1000.
		    path_in_he = 0.
		    path_in_he_after = 0.
		    rear_wall_zpos = -1.5 + tar_zoff(1) 
		    front_wall_zpos = 1.5 + tar_zoff(1) 
		  if ( z .lt. rear_wall_zpos) then ! before target cell rear wall
		    if ( z .ge. (rear_wall_zpos-.5) ) then ! In He bfore cell
		       path_in_he = (rear_wall_zpos - z)/cos(th_spec-dydz)
		    else
		       path_in_he = .5/cos(th_spec-dydz)
		    endif
                    if ( path_in_he .gt. 0) then
                      call enerloss_new(path_in_he,0.145d0,
     >2.0d0,4.0d0,ep,.511d0,1,temp_eloss) ! He before target
		    tot_ran_escat_eloss  = tot_ran_escat_eloss + temp_eloss
		    endif
		    x_rear_wall = x + (rear_wall_zpos-z)*tan(th_spec-dydz)
		     if ( x_rear_wall .lt. 1.27) then ! hit cell		       
		       x_front_wall = x + (front_wall_zpos-z)*tan(th_spec-dydz)
		       if ( x_front_wall .lt. 1.27) then 
                          path = sqrt( (x_front_wall-x_rear_wall)**2+ 3.0*3.0)
		          path_in_he_after = .5/cos(th_spec-dydz)
                       else
			  z_side_wall = rear_wall_zpos+(1.27-x_rear_wall)/tan(th_spec-dydz) 
                          path = sqrt( (1.27-x_rear_wall)**2+ z_side_wall**2) 
			  path_in_he_after = .5/cos(th_spec-dydz) + 
     >  sqrt( (x_front_wall-1.27)**2+ (front_wall_zpos-z_side_wall)**2)
		       endif
		       path=path/2.
                     else
			path_in_he_after = 3.5/cos(th_spec-dydz)
		     endif
                  elseif ( z .ge. rear_wall_zpos .and. z .le. front_wall_zpos) then
		       x_rear_wall = -1000.
		       x_front_wall = x + (front_wall_zpos-z)*tan(th_spec-dydz)
		       if ( x_front_wall .lt. 1.27) then 
                          path = sqrt( (x_front_wall-x)**2+ (front_wall_zpos-z)**2)
                          path_in_he_after = .5/cos(th_spec-dydz)
                       else
			  z_side_wall =(1.27-x)/tan(th_spec-dydz)
                          path = sqrt( (1.27-x)**2+ z_side_wall**2) 
			  path_in_he_after = .5/cos(th_spec-dydz) + 
     >  sqrt( (x_front_wall-1.27)**2+ (front_wall_zpos-z_side_wall)**2)
		       endif
		       path=path/2. ! assum PF 50% split path between (NH or ND) and He
		  elseif ( z .gt. front_wall_zpos .and. z .le. 2.0) then
		       path_in_he_after = (front_wall_zpos + .5 - z)/cos(th_spec-dydz)		       
                  endif                  
		  if ( path .ne. 0) then
                       call enerloss_new(path,tar_den(tartype),
     >tar_z(tartype),tar_a(tartype),ep,.511d0,1,temp_eloss) ! NH or Nd in target assume 50% PF
		       tot_ran_escat_eloss  = tot_ran_escat_eloss + temp_eloss
                       call enerloss_new(path,0.145d0,
     >2.0d0,4.0d0,ep,.511d0,1,temp_eloss) ! He in target assume 50% PF
		       tot_ran_escat_eloss  = tot_ran_escat_eloss + temp_eloss
		       if ( path_in_he_after .ne. 0) then
                       call enerloss_new(path_in_he_after,0.145d0,
     >2.0d0,4.0d0,ep,.511d0,1,temp_eloss) ! 
		       tot_ran_escat_eloss  = tot_ran_escat_eloss + temp_eloss
		       endif
		  endif
		  path_tail = .01016/cos(th_spec-dydz)
                  call enerloss_new(path_tail,2.7d0,
     >13.0d0,27.0d0,ep,.511d0,1,temp_eloss) ! tailpiece 
		  tot_ran_escat_eloss  = tot_ran_escat_eloss + temp_eloss
		  path_4k = .00254/cos(th_spec-dydz)
                  call enerloss_new(path_4k,2.7d0,
     >13.0d0,27.0d0,ep,.511d0,1,temp_eloss) ! 4k
		  tot_ran_escat_eloss  = tot_ran_escat_eloss + temp_eloss
		    mscat_len = (path_in_he+path_in_he_after+path)*0.145/94.3 
     >                 + path*.917/43.255 + (path_tail + path_4k)*2.7/24.0
               endif
         mscat_len = mscat_len + .040*2.7/24.01 + 15*.00121/36.66
     >         + .015*2.54*.74/55.2 + .005*2.54*1.39/39.95
                  call enerloss_new(.040d0,2.7d0,
     >13.0d0,27.0d0,ep,.511d0,1,temp_eloss) ! scattering chamber window
		  tot_ran_escat_eloss  = tot_ran_escat_eloss + temp_eloss
                  call enerloss_new(15.d0,0.00121d0,
     >7.2d0,14.4d0,ep,.511d0,1,temp_eloss) ! air
		  tot_ran_escat_eloss  = tot_ran_escat_eloss + temp_eloss
                  call enerloss_new(0.015d0*2.54d0,0.74d0,
     >2.67d0,4.67d0,ep,.511d0,1,temp_eloss) ! kevlar
		  tot_ran_escat_eloss  = tot_ran_escat_eloss + temp_eloss
                  call enerloss_new(0.005d0*2.54d0,1.39d0,
     >4.545d0,8.735d0,ep,.511d0,1,temp_eloss) ! mylar
		  tot_ran_escat_eloss  = tot_ran_escat_eloss + temp_eloss

c       
	  ep=ep - tot_ran_escat_eloss 
c
	  enu  = ebeam_ran - ep
	  ebeam = ebeam_ran
	  
C radiative corrections

          rad_weight = 1.0
          erad  = 0.0
          eprad = 0.0
          success = .false.
	    etemp=ebeam
	    eptemp=ep
          if(radiate) then
	    call radc_init(zz(imat))
c	    write(*,*) trial,' call radc_init_ev'
            call radc_init_ev(tb,ta,etemp,eptemp,th_ev,ue)
            call generate_rad(elastic_event,etemp,eptemp
     >     ,ep_rad_max,ep_rad_min,ue,rad_weight,success,erad,eprad)
          endif

C calculate kinematics  
c calculate vertex quantities for QFS xn

	  ebeam_vertex = ebeam - erad 
	  eprime_vertex = ep
c 
 	  dpps = dpps - eprad/p_spec*100.
	  ep = ep -eprad
	  if ( elastic_event .eq. 1 .and. erad .gt. 0 ) then
	     eprime_vertex = mp*ebeam_vertex/(ebeam_vertex*(1-ue(3))+mp)
	     dpps = (eprime_vertex-p_spec)/p_spec*100.
	     ep = eprime_vertex
	  endif   
	  enu_vertex = ebeam_vertex -eprime_vertex
          if(enu_vertex .le. 0) go to 500
 	  Q2_vertex  = 4.0*1d-6*ebeam_vertex*eprime_vertex*sin(th_ev/2)**2
          m_n = Mn*1.0d-3
          W_vertex   = 2.*m_n*enu_vertex*1.d-3 + m_n**2 - Q2_vertex
	  if ( W_vertex .le. 0) then
	     go to 500
          else
	     W_vertex = sqrt(W_vertex)
	  endif
	  if ( a(imat) .eq. 1 .and. zz(imat) .eq. 1 
     >       .and. elastic_event .eq. 0) then ! inelastic ep event
	     if ( W_vertex .le. 1.073) go to 500
          endif
c
c
C save init values for later

	  xbeam_init = x
	  ybeam_init = y
	  zbeam_init = z
	  xtar_init = xs
	  ytar_init = ys
	  ztar_init = zs
	  dpp_init = dpps
	  dth_init = dydzs
	  dph_init = dxdzs

C transport through spectrometer and do reconstruction if asked
	  dpp_recon  = 0.
	  dth_recon  = 0.	
          dph_recon  = 0.	
	  xtar_recon = 0.
	  ytar_recon = 0.
	  ztar_recon = 0.

          ok_spec = .true.
c  multiple scattering
	  p_temp = (dpps/100.+1)*p_spec
	  if ( mscat_len .ne. 0) then
           call musc(me2,p_temp,mscat_len,dxdzs,dydzs) 
c	  write(*,*) ' mscat imat = ',imat,z,mscat_len,ta,dxdzs,dph_init,dydzs,dth_init
	  endif
c	  
	  if ( fieldon) then
c
c code tracks particle to z=100 
c
	  call track_from_tgt(xs, ys, zs, dxdzs, dydzs, -ep,me,-1,ok_spec)
          if (.not. ok_spec) nfail_track_from_tgt = nfail_track_from_tgt + 1
c
C drift back to zs = 0, the plane through the target center
c
	  xs = xs - zs * dxdzs
	  ys = ys - zs * dydzs
	  zs = 0.0d00
c
	  endif
c
          if (ok_spec) then
	  call mc_hms(p_spec, th_spec, dpps, xs, ys, zs, dxdzs, dydzs, 
     >                x_fp, dx_fp, y_fp, dy_fp, m2, ms_flag, wcs_flag, 
     >                decay_flag, resmult, y, x, ok_spec)
          endif
          if(ok_spec) then
  	    dpp_recon  = dpps
	    dth_recon  = dydzs
            dph_recon  = dxdzs
	    xtar_recon = xs
	    ytar_recon = ys
c mkj use xbeam_init
            ztar_recon = (ytar_recon*cos_ts+xbeam_init)
     >                /tan(th_spec-dth_recon)+ytar_recon*sin_ts
          endif
c
	  cos_ev_recon = (cos_ts+dydzs*sin_ts)/sqrt(1+dydzs**2+dxdzs**2)
	  th_ev_recon  = acos(cos_ev_recon)                            
	  theta_ev_recon = acos((cos_ts+dydzs*sin_ts)/sqrt(1+dydzs**2)) 
          phi_ev_recon = asin(dxdzs/sqrt(1+dxdzs**2+dydzs**2))   
	  ep_recon = p_spec*(1+0.01*dpp_recon)
	  hsphi = atan2((dydzs*cos_ts-sin_ts),dxdzs)
          hsphi  =  hsphi + 270.*3.14159/180.
c
          beta_temp  = 1./sqrt(1.+(.511d0/ep_recon)**2)
          gamma_temp = sqrt(1.+(.511d0/ep_recon)**2)/(.511d0/ep_recon)
          velocity   = log(beta_temp*gamma_temp)/log(10.)
	  tot_mp_escat_eloss = 0
	  if (tartype .eq. 5 .or. tartype .eq. 6) then
	  th_temp=tar_len(tartype)*tar_den(tartype)/2./cos(th_ev_recon)
	  call loss(0,tar_z(tartype),tar_a(tartype),th_temp,tar_den(tartype),velocity,temp_eloss) ! back_loss
	  tot_mp_escat_eloss = tot_mp_escat_eloss  + temp_eloss
	  call loss(0,2.d0,4.d0,0.184d0,0.145d0,velocity,temp_eloss) ! back_he_loss
	  tot_mp_escat_eloss = tot_mp_escat_eloss + temp_eloss
	  call loss(0,13.d0,27.d0,0.01056d0,2.7d0,velocity,temp_eloss) ! cell_wall_loss
	  tot_mp_escat_eloss = tot_mp_escat_eloss + temp_eloss
	  elseif ( tartype .eq. 1) then
	     th_temp = .145*abs(4-tar_len(tartype))/2./cos(th_ev_recon)
	     call loss(0,2.d0,4.d0,th_temp,0.145d0,velocity,temp_eloss) 
	     tot_mp_escat_eloss =temp_eloss
	     th_temp=tar_len(tartype)*tar_den(tartype)/2.
	     call loss(0,tar_z(tartype),tar_a(tartype),th_temp,tar_den(tartype),velocity,temp_eloss) 
	     tot_mp_escat_eloss = tot_mp_escat_eloss + temp_eloss
	  else
	     th_temp=tar_len(tartype)*tar_den(tartype)/2.
	     call loss(0,tar_z(tartype),tar_a(tartype),th_temp,tar_den(tartype),velocity,temp_eloss) 
	     tot_mp_escat_eloss =temp_eloss
	  endif

c
	  call loss(0,13.d0,27.d0,0.08915d0,2.7d0,velocity,temp_eloss) ! scat_win_loss
	  tot_mp_escat_eloss = tot_mp_escat_eloss + temp_eloss
	  call loss(0,7.32d0,14.68d0,0.018d0,0.00121d0,velocity,temp_eloss) ! air_loss
	  tot_mp_escat_eloss = tot_mp_escat_eloss + temp_eloss
	  call loss(0,2.67d0,4.67d0,0.0491d0,0.87864d0,velocity,temp_eloss) ! h_win_loss
	  tot_mp_escat_eloss = tot_mp_escat_eloss + temp_eloss
	  ep_recon = ep_recon + tot_mp_escat_eloss  
c
 	  Q2_recon  = 4.0*1d-6*ebeam_mp*ep_recon*sin(th_ev_recon/2)**2
          m_n = Mn*1.0d-3
          W_recon   = 2.*m_n*(ebeam_mp-ep_recon)*1.d-3 + m_n**2 - Q2_recon
	  if ( W_recon .le. 0) then
	     W_recon = 0
	  else
	     W_recon = sqrt(W_recon)
	  endif
C calculate cross sections and yield normalization factor

	  th = th_ev*180/acos(-1.)
c
          call qfs(zz(imat), a(imat), ebeam_vertex, th, enu_vertex, delnu, 
     >             rad_qfs,sigma_qfs, elastic_event,scalein)
	  if ( elastic_event .eq. 1) then
	       sigrad = sigma_qfs*1.0d33  !convert cm**2/sr to nb/sr
          else
	       sigrad = sigma_qfs*1.0d33  !convert cm**2/MeV/sr to nb/MeV/sr
	  endif
          sig_pb = 0
	  if ( a(imat) .gt. 2) then
            call pb_ext_sub(zz(imat),a(imat),ebeam_vertex,enu_vertex,th,sig_pb)
	    sig_pb = sig_pb/1000.
          endif
	  if(ok_spec) then
	     if ( a(imat) .eq. 1 .and. elastic_event .eq. 1) then
	    normrate = sigrad*rad_weight*proton_wgt*
     > lumin(imat)*domega/ver(imat)/n_trials
	    else
	    normrate = sigrad*rad_weight*proton_wgt*
     > lumin(imat)*denergy*domega/ver(imat)/n_trials
            endif	       
            normrate2 = normrate
            if ( a(imat) .gt. 2) then 
	    normrate2 = sig_pb*rad_weight*proton_wgt*
     > lumin(imat)*denergy*domega/ver(imat)/n_trials
            endif
	  else
	    normrate=0.0
            normrate2 = 0.0
	  end if

c mkj temp comment out.
c 
	  cer_eff=1.
	  dc_eff=1.
c	   call cer_effcorr(dpp_recon,cer_eff)
c	   call dc_effcorr(x_fp,dc_eff)


C write ntuple if particle made it

	  if(ok_spec) then
	    hut(1)  = normrate
	    hut(2)  = a(imat)
	    hut(3)  = W_recon
	    hut(4)  = dpp_init
	    hut(5)  = dpp_recon
	    hut(6)  = th_ev_init*degrad
	    hut(7)  = th_ev_recon*degrad
	    hut(8)  = theta_ev_init*degrad
	    hut(9)  = theta_ev_recon*degrad
	    hut(10) = hsphi_init*degrad
	    hut(11) = hsphi*degrad
	    hut(12) = dph_init
	    hut(13) = dph_recon
	    hut(14) = dth_init
	    hut(15) = dth_recon
	    hut(16) = Q2_recon
	    hut(17) = x_fp
	    hut(18) = y_fp
	    hut(19) = dx_fp
	    hut(20) = dy_fp
	    hut(21) = xtar_init
	    hut(22) = xtar_recon
	    hut(23) = ytar_init
	    hut(24) = ytar_recon
	    hut(25) = ztar_recon
	    hut(26) = zbeam_init
	    hut(27) = ybeam_init
	    hut(28) = xbeam_init
	    hut(29) = elastic_event
	    call hfn(1,hut)
	  end if

	  ave_rate = 0
	  do i=1, nmat
	     if ( trials(i) .gt. 0) then
	    ave_rate = ave_rate+rate(i)/trials(i)
             endif
	  end do

	  trials(imat) = trials(imat) + 1
	  if (ok_spec) then
	    sigma_hms(imat) = sigma_hms(imat) + sigrad*rad_weight
	    if (a(imat) .eq. 1 .and. elastic_event .eq. 1) then
	    rate(imat)=rate(imat)+sigrad*lumin(imat)*domega*rad_weight
	    else
	    rate(imat)=rate(imat)+sigrad*lumin(imat)*domega*denergy*rad_weight
	    endif
c	    write(*,*) rate(imat),sigrad,lumin(imat),rad_weight
	  end if

C cut on reconstructed quantities

	  if ((abs(dpp_recon).gt.cut_dpp) .or.
     >	      (abs(dth_recon).gt.cut_dth) .or.
     >	      (abs(dph_recon).gt.cut_dph) .or.
     >	      (abs(ztar_recon).gt.cut_z)) then
	    goto 500		!finish loop if failed
	  endif

C compute sums for calculating reconstruction variances

	  dpp_var(1) = dpp_var(1) + (dpp_recon - dpp_init)
	  dth_var(1) = dth_var(1) + (dth_recon - dth_init)
	  dph_var(1) = dph_var(1) + (dph_recon - dph_init)
	  ztg_var(1) = ztg_var(1) + (ztar_recon - ztar_init)

	  dpp_var(2) = dpp_var(2) + (dpp_recon - dpp_init)**2
	  dth_var(2) = dth_var(2) + (dth_recon - dth_init)**2
	  dph_var(2) = dph_var(2) + (dph_recon - dph_init)**2
	  ztg_var(2) = ztg_var(2) + (ztar_recon - ztar_init)**2

C we are done with this event, whether GOOD or BAD

          if(ok_spec) then
            ntp_events = ntp_events + 1
            if ( ntp_events .eq. ntp_events_max ) then
cc                                  [close the current ntuple and open new one]
               call hrout(1,i,' ')
	       call hrend('HUT')                  
               print *, 'CLOSE the current ntuple:', hbname

	       ntp_events = 0  !! initialize ntuple counter 
	       k=lnblnk(hbname0)
	       ntp_number = ntp_number + 1     !! Next ntuple file number
	       call DIGIT2STRING(ntp_number, ntp_number_char)
	       l=lnblnk(ntp_number_char)
	       hbname=hbname0(1:k)//'.'//ntp_number_char(1:l)//'.hbook'
               print *, 'NEW ntuple file =', hbname
	       call hropen(30,'HUT',hbname,'NQ',1024,i)
	       call hbookn(1,'HUT NTUPLE',nvars,'HUT',1024,hut_nt_names)
	    endif 

	  endif

500	  continue
	enddo				!End of M.C. loop

C----------------------------------------------------------------------------C
C                         End of Monte-Carlo loop                            C
C----------------------------------------------------------------------------C

C close NTUPLE file

 533	call hrout(1,i,' ')
	call hrend('HUT')

	write (chanout,1002)
	write (chanout,1003) e,p_spec,th_spec*degrad
        write (chanout,1004) dpp_lo,dpp_hi,th_lo*1000.,th_hi*1000.,ph_lo
     > *1000.,ph_hi*1000.,raster_radius,raster_xoff,raster_yoff,packfrac
	write (chanout,2004) p_flag,ms_flag,wcs_flag,rad_qfs,fieldon,th_b,recon,radiate
	write (chanout,1005) n_trials

C indicate where particles are lost in spectrometer

	write (chanout,1015)
     >	nfail_track_from_tgt,hSTOP_slit_hor,hSTOP_slit_vert,hSTOP_slit_oct,
     >	hSTOP_Q1_in,hSTOP_Q1_mid,hSTOP_Q1_out,
     >	hSTOP_Q2_in,hSTOP_Q2_mid,hSTOP_Q2_out,
     >	hSTOP_Q3_in,hSTOP_Q3_mid,hSTOP_Q3_out,
     >	hSTOP_D1_in,hSTOP_D1_out,nfail_track_to_tgt

	write (chanout,1006)
     >	hSTOP_trials,hSTOP_hut,hSTOP_dc1,hSTOP_dc2,hSTOP_scin,hSTOP_cal,
     >  hSTOP_successes,hSTOP_successes

C compute reconstruction resolutions

	write(6,*) hSTOP_trials,' Trials',hSTOP_successes,' Successes'

	if (hSTOP_successes.eq.0) hSTOP_successes=1 !cheat to avoid / by 0
	t1=sqrt(max(0.,dpp_var(2)/hSTOP_successes-
     >              (dpp_var(1)/hSTOP_successes)**2))
	t2=sqrt(max(0.,dth_var(2)/hSTOP_successes-
     >              (dth_var(1)/hSTOP_successes)**2))
	t3=sqrt(max(0.,dph_var(2)/hSTOP_successes-
     >              (dph_var(1)/hSTOP_successes)**2))
	t4=sqrt(max(0.,ztg_var(2)/hSTOP_successes-
     >              (ztg_var(1)/hSTOP_successes)**2))

	write (chanout,1011) dpp_var(1)/hSTOP_successes,t1,
     >        dth_var(1)/hSTOP_successes,t2,dph_var(1)/hSTOP_successes,
     >        t3,ztg_var(1)/hSTOP_successes,t4

	write(6,1011) dpp_var(1)/hSTOP_successes,t1,dth_var(1)/hSTOP_successes,
     >	      t2,dph_var(1)/hSTOP_successes,t3,ztg_var(1)/hSTOP_successes,t4

C we are done!	

	rate_hms = 0
	write(chanout,*) '# trials and rate for each target layer:'
	print *,'Target layer, num trials, rate/trials'
	do i=1, nmat
	   write(chanout,1020) trials(i), rate(i)/trials(i)
	   print '(i3,1x,f10.0,1x,g18.8)',i, trials(i), rate(i)/trials(i)
	   rate_hms = rate_hms + rate(i)/trials(i)
	end do

	write(chanout,*) ' '
	write(chanout,*) 'Total rate: ', rate_hms

	print *, 'estimated rate:', rate_hms          ! *domega*denergy
	print *, 'domega=', domega
	print *, 'de = ', denergy
	print *, 'de*denergy=', domega*denergy

	stop

C---------------------------- Format Statements -----------------------------

1001	format(a)
1002	format('!',/,'! Uniform illumination Monte-Carlo results')
1003	format(g18.8,' =  Beam Energy (MeV)',/,
     >  '!',/'! Spectrometer setting:',/,'!',/,
     >	g18.8,' =  P  spect (MeV)',/,
     >	g18.8,' =  TH spect (deg)')

1004	format('!',/'! Monte-Carlo limits:',/,'!',/,
     >	g18.8,' =  dP/P Lower Limit                                 (%)',/,
     >	g18.8,' =  dP/P Upper Limit                                 (%)',/,
     >	g18.8,' =  Theta_scatter, Lower Limit                      (mr)',/,
     >	g18.8,' =  Theta_scatter, Upper Limit                      (mr)',/,
     >	g18.8,' =  Phi_scatter, Lower Limit                        (mr)',/,
     >	g18.8,' =  Phi_scatter, Upper Limit                        (mr)',/,
     >	g18.8,' =  Raster Radius                                   (cm)',/,
     >	g18.8,' =  Horizontal raster offset            (+=beam left,cm)',/,
     >	g18.8,' =  Vertical raster offset                     (+=up,cm)',/,
     >	g18.8,' =  packing fraction                    ')

2004	format('!',/'! Flags:',/,'!',/,
     >	i1,' =   particle identification: e=0, p=1, d=2, pi=3, ka=4 ',/,
     >	l1,' =   flag for multiple scattering                       ',/,
     >	l1,' =  flag for wire chamber smearing                      ',/,
     >	l1,' =  flag for QFS internal radiation: 1=on, 0=off     ',/,
     >	l1,' =  flag to turn on target magnetic field: 1=on, 0=off',/,
     >	g18.8,' =  target magnetic field direction: minus sign is beam left',/,
     >	l1,' =  reconstruction flag: 1=yes, 0=no                  ',/,
     >	l1,' =  radiation flag (simc way): 1=on, 0=off')

1005	format('!',/,'! Summary:',/,'!',/,i12,' Monte-Carlo trials:')

1006	format(i12,' Initial Trials',/
     >	i12,' Trials cut in the hut',/
     >	i12,' Trials cut in dc1',/
     >	i12,' Trials cut in dc2',/
     >	i12,' Trials cut in scin',/
     >	i12,' Trials cut in cal',/
     >	i12,' Trials made it thru the detectors and were reconstructed',/
     >	i12,' Trials passed all cuts and were histogrammed.',/
     >	)

1011	format(
     >	'DPP ave error, resolution = ',2g18.8,' %',/,
     >	'DTH ave error, resolution = ',2g18.8,' mr',/,
     >	'DPH ave error, resolution = ',2g18.8,' mr',/,
     >	'ZTG ave error, resolution = ',2g18.8,' cm')


1015	format(/,
     >	i12,' failed track from target',/
     >	i12,' stopped in the FIXED SLIT HOR',/
     >	i12,' stopped in the FIXED SLIT VERT',/
     >	i12,' stopped in the FIXED SLIT OCTAGON',/
     >	i12,' stopped in Q1 ENTRANCE',/
     >	i12,' stopped in Q1 MIDPLANE',/
     >	i12,' stopped in Q1 EXIT',/
     >	i12,' stopped in Q2 ENTRANCE',/
     >	i12,' stopped in Q2 MIDPLANE',/
     >	i12,' stopped in Q2 EXIT',/
     >	i12,' stopped in Q3 ENTRANCE',/
     >	i12,' stopped in Q3 MIDPLANE',/
     >	i12,' stopped in Q3 EXIT',/
     >	i12,' stopped in D1 ENTRANCE',/
     >	i12,' stopped in D1 EXIT',/
     >	i12,' failed track to target',/
     >	)

 1020	format(f10.0,1x,g18.8)


 9999	print *, 'Usage: mc_hms_single [input filename] '
     > //'[target filename] [kinematics specs]'
	print *, 'Input files are located in the infiles/ directory.'
	print *, 'Target files are located in the tarfiles/ directory.'
	stop
	end


      SUBROUTINE DIGIT2STRING(INP,  FILN)
C*************************************************************      
C
C     ... Produce character values for each digit.
C
C     INPUT : INP     - the number of ntuple file
C     
C     OUTPUT: FILN    - string with the ntuple number
C     
C*************************************************************      
C     
      INTEGER   LENGTH, I,J,K,INP
      CHARACTER FILN*80, DIGIT*1
C           
      I = IABS( INP )
      LENGTH = 0
      FILN = ' '
    1 LENGTH = LENGTH + 1
        J = I/10
        K = I - 10*J
        DIGIT  = CHAR( 48+K )
        FILN = DIGIT // FILN
        I = J
      IF (I .GT. 0) GO TO 1
      
      RETURN
      END
