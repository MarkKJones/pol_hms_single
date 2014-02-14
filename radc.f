!--------------------------------------------------------------------------
	subroutine generate_rad(elastic_event,Ein,Ep,E_max,E_min,ue,
     >                          rad_weight,success,E_rad,Ep_rad)

	implicit none
c	include 'simulate.inc'
	include 'radc.inc'
	include 'constants.inc'
	
	integer i, peaked_basis_flag
	real*8 x, rad_weight, ue(3)
	real*8 peaked_rad_weight, basicrad_val_reciprocal
	real*8 basicrad_weight, extrad_phi
        real*8 Ein, Ep, E_max, E_min, E_rad, Ep_rad
	real*8 elastic_event
	logical success

	real*8 grnd
	real*8 ebeam_min,ebeam_max
	
!---Generate Radiation-----------------------------------------------------

! Initialize.

	success = .false.
	rad_weight = 1
	do i = 1,2
	  Egamma_used(i) = 0.0
	enddo

! Which tail has to be radiated? Set ntail=0 if all tails allowed to radiate.

	peaked_basis_flag = 1
	if (rad_flag.le.1) then
	  peaked_basis_flag = 0
	  x = grnd()
	  if (x.ge.frac(1)) then
	    ntail = 2
	  else
	    ntail = 1
	  endif
	else if (rad_flag.eq.2) then
	  ntail = int(grnd()*2.) + 1
	  if (ntail.eq.3) ntail=2
	else if (rad_flag.eq.3) then
	  ntail = 0
	else
	  stop 'Idiot! rad_flag is set stupidly'
	endif

! RADIATE TAIL #1: the incident electron tail

	if (doing_tail(1) .and. (ntail.eq.0.or.ntail.eq.1)) then
	  Egamma_min(1) = 0.
	  Egamma_max(1) = Ein-E_min
	  if ( elastic_event .ge. 1) then
	     ebeam_max = mp*e_max/(mp-e_max*(1-ue(3)))
	     ebeam_min = mp*e_min/(mp-e_min*(1-ue(3)))
	     egamma_min(1) = Ein - ebeam_max
	     egamma_max(1)  = Ein - ebeam_min
	     egamma_max(1) =min(egamma_max(1),e_max)
	  endif

! ... radiate  

	  call basicrad(1*peaked_basis_flag,Egamma_min(1),Egamma_max(1),
     >          Egamma_used(1),basicrad_weight,basicrad_val_reciprocal)
	  if (basicrad_weight.le.0) then
	     return
	  endif
          Ein = Ein - Egamma_used(1)
	  if ( elastic_event .eq. 1) Ep = mp*Ein/(Ein*(1-ue(3))+mp)
	  rad_weight = rad_weight * basicrad_weight
	  if (rad_flag.le.1) then
	    rad_weight = peaked_rad_weight(Ep,ue,Ein,Egamma_used(1),Egamma_min(1),
     >		Egamma_max(1),basicrad_val_reciprocal,basicrad_weight)
	  else
	    rad_weight=rad_weight*extrad_phi(1,Ein,Ep,Egamma_used(1))
	  endif
	endif		!end of tail 1 (incoming electron)

! RADIATE TAIL #2: the scattered electron tail

	if (doing_tail(2) .and. (ntail.eq.0.or.ntail.eq.2)) then
	  Egamma_min(2) = 0.
	  Egamma_max(2) = Ep - E_min
! ... radiate
	  call basicrad(2*peaked_basis_flag,Egamma_min(2),Egamma_max(2),
     >          Egamma_used(2),basicrad_weight,basicrad_val_reciprocal)
	  if (basicrad_weight.le.0) then
	    return
	  endif
	  rad_weight = rad_weight * basicrad_weight
	  if (rad_flag.le.1) then
	    rad_weight = peaked_rad_weight(Ep,ue,Ein,Egamma_used(2),Egamma_min(2),
     >		Egamma_max(2),basicrad_val_reciprocal,basicrad_weight)
	  else 
	    rad_weight=rad_weight*extrad_phi(2,Ein,Ep,Egamma_used(2))
	  endif
	endif		!end of tail 2 (outgoing electron)

! ... remove radiation from incoming electron.

	E_rad  = Egamma_used(1)
        Ep_rad = Egamma_used(2)

! Complete determination of the event weight due to radiation.

	rad_weight = rad_weight/hardcorfac
	success = .true.

	return
	end

!--------------------------------------------------------------------------
	subroutine basicrad(itail,Egamma_lo,Egamma_hi,Egamma,
     >                      weight,val_reciprocal)

	implicit none
	include 'radc.inc'

	integer itail
	real*8 Egamma_lo, Egamma_hi, Egamma, weight, val_reciprocal
	real*8 power_lo, power_hi
	real*8 x, y, ymin

	real*8 grnd

!--------------------------------------------------------------------------
!
! Generate Egamma inside the requested region, using the function 
!	    g * Egamma**(g-1)
!	-------------------------------
!	(Egamma_max**g - Egamma_min**g)
! which has the BASICRAD shape, is "invertible" as a probability 
! distribution,
! and integrates to 1 over the requested region { Egamma_min to Egamma_max }.
! Each radiative event we generate with this then has to be weighted with
! the _actual_ radiative tail form (which doesn't necessarily integrate to 1)
! _divided_by_ this generating function, so we return the RECIPROCAL of the 
! function's value (at the selected Egamma) as VAL_RECIPROCAL. If we're not 
! using anything fancier than BASICRAD, we can immediately determine what
! that 
! quotient is:
!		  C
!	weight = --- * (Egamma_max**g - Egamma_min**g)
!		  g
! We return that as (you guessed it!) WEIGHT.
! If we want to use BASICRAD x the phi function (correction to external tail
! shape), all we need to do is multiply WEIGHT by PHI evaluated at the
! selected photon energy. That's dont outside this routine, in GENERATE_RAD.
! Note that if we're not using the BASICRAD prescription at all (except as
! a generation tool), WEIGHT is not used again. 
!
! Note this routine handles the possibilities Egamma_lo < 0, Egamma_hi <
! Egamma_lo, and Egamma_hi < 0
!---------------------------------------------------------------------------

! Initialize

	Egamma = 0.0
	weight = 0.0
	val_reciprocal = 0.0

! ... is radiation in this direction turned off?

	if(itail.eq.0) itail=4
	if (g(itail).le.0) then
	  weight = 1.0
	  return
	endif

! ... do we have a chance?

	if (Egamma_hi.le.Egamma_lo .or. Egamma_hi.le.0) return

! ... want to evaluate powers as INfrequently as possible!

	power_hi = Egamma_hi**g(itail)
	power_lo = 0.0
	if (Egamma_lo.gt.0) power_lo = Egamma_lo**g(itail)

! ... obtain y region: from (Egamma_lo/Egamma_hi)**g to 1, generate flat y.

	ymin = power_lo/power_hi
	y = ymin + grnd()*(1.-ymin)
	x = y**(1./g(itail))
	Egamma = x*Egamma_hi

! The value of our generating function at the selected point

	if (Egamma.gt.0) val_reciprocal = Egamma**(1.-g(itail)) *
     >		(power_hi-power_lo) / g(itail)

! Event weight = probability radiation falls inside requested range,
! in full BASICRAD prescription

	weight = c(itail)/g(itail) * (power_hi-power_lo)

	if(itail.eq.4) itail=0

	return
	end

!--------------------------------------------------------------------------

	real*8 function peaked_rad_weight(Ep,ue,Ein,Egamma,
     >		emin,emax,basicrad_val_reciprocal,basicrad_weight)

	implicit none
	include 'radc.inc'
c	include 'structures.inc'

	real*8 Egamma, basicrad_val_reciprocal, basicrad_weight
	real*8 t1, t2, r, phi_ext, ein, eout, emin, emax
	real*8 bremos, extrad_phi, gamma !,brem
	real*8 dhard, dsoft_intmin, dsoft_intmax
	real*8 dsoft_ext1, dsoft_ext2, dsoft_extmin
	real*8 dsoft_int_primemin, dsoft_int_primemax
	real*8 dsoft_ext1_prime, dsoft_ext2_prime
	real*8 dsoft_ext_prime, dsoft_extmax, eul
	real*8 dsoft_int,dsoft_ext,dsoft_int_prime
        real*8 Ep, ue(3)

	real*8 zero
	parameter (zero=0.0d0)	!double precision zero for subroutine calls

	basicrad_val_reciprocal=basicrad_val_reciprocal+0. !avoid unused variable error
! Compute a more precise value for the radiative probability at the
! selected Egamma, more precise than the BASICRAD distribution used to 
! generate Egamma that is.
!
! NB: This subroutine only deals with calculations done in the PEAKING
! APPROXIMATION basis -
!    i.e. we only get here if RAD_FLAG = 0 or 1
!
! ... (sort of) KLUGE -- if Egamma < res limit, these things might blow up,
! ... so just screw the correction and leave
	ein = Ein
	eout= Ep
	eul = 0.577215665

! ... External

	phi_ext = 1.0
	if (extrad_flag.le.2) then
	  phi_ext = extrad_phi(0,ein,eout,Egamma)
	  if (rad_flag.eq.1) then
	    peaked_rad_weight = basicrad_weight * phi_ext
	    return
	  endif
	  if(emin.gt.0)then
	    dsoft_extmin = log(g_ext/c_ext(0)/emin**g_ext)
	  endif
	  dsoft_extmax = log(g_ext/c_ext(0)/emax**g_ext)
	  dsoft_ext_prime = -g_ext/Egamma
	else
	  t1 = bt(1)/etatzai
	  t2 = bt(2)/etatzai
	  call extrad_friedrich(ein,Egamma,t1,dsoft_ext1,dsoft_ext1_prime)
	  call extrad_friedrich(eout,Egamma,t2,dsoft_ext2,dsoft_ext2_prime)
	  dsoft_ext = dsoft_ext1 + dsoft_ext2
	  dsoft_ext_prime = dsoft_ext1_prime + dsoft_ext2_prime
	endif

! ... Internal
! ........ use full BREM calculation of deltas

	if (rad_flag.eq.0) then
	    if(emin.gt.0)then
	      r = bremos(emin, zero, zero, ein, Ep*ue(1), Ep*ue(2),Ep*ue(3),
     >			dsoft_intmin, dhard, dsoft_int_primemin)
	    else 
	      dsoft_intmin=1.0
	    endif
	      r = bremos(emax, zero, zero, ein, Ep*ue(1), Ep*ue(2),Ep*ue(3),
     >			dsoft_intmax, dhard, dsoft_int_primemax)

! ........ use basic calculation of internal deltas

	else
	  dsoft_int = log(g_int/c_int(0)/Egamma**g_int)
	  dsoft_int_prime = -g_int/Egamma
	endif

! ... All together now

	if(emin.gt.0)then
	  peaked_rad_weight = c_ext(0)/g_ext*(exp(-dsoft_intmax)*
     >		emax**g_ext-exp(-dsoft_intmin)*emin**g_ext)
	else
	  peaked_rad_weight = c_ext(0)/g_ext*(exp(-dsoft_intmax)*emax**g_ext)
	endif
	peaked_rad_weight = peaked_rad_weight * exp(-eul*g(4))/gamma(1.+g(4))
     >		* gamma(1.+g(4)-bt(1)-bt(2))*gamma(1.+bt(1))
     >		* gamma(1.+bt(2))/gamma(1.+g(4))
	if (peaked_rad_weight.lt.0) peaked_rad_weight = 0

	return
	end

!---------------------------------------------------------------------------

	subroutine extrad_friedrich(Ei,Ecutoff,trad,dbrem,dbrem_prime)

	implicit none
c	include 'simulate.inc'
	include 'radc.inc'

	real*8	Ei, Ecutoff, trad, x, dbrem, dbrem_prime

	x = Ecutoff/Ei
	dbrem = trad*(-(etatzai-0.5) - etatzai*log(x) + etatzai*x - 0.5*x**2)
	dbrem_prime = -trad/Ei * (etatzai/x - etatzai + x)

	return
	end

!--------------------------------------------------------------------------

	real*8 function extrad_phi(itail,E1,E2,Egamma)

	implicit none
	include 'radc.inc'

	integer itail
	real*8 E1, E2, Egamma, E(2), x, t, gamma

	E(1) = E1
	E(2) = E2

! Compute the multiplicative external correction function PHI

! ... CASE 1: phi = 1

	extrad_phi = 1.0

! ... CASE 2: phi correctly computed according to Dave's formulas

	if (extrad_flag.eq.2) then
	  if (itail.eq.0) then
	    extrad_phi = 1. - (bt(1)/E(1)+bt(2)/E(2)) / (g(1)+g(2)) * Egamma
	  else if (itail.eq.1.or.itail.eq.2) then
	    extrad_phi = 1. - bt(itail)/E(itail) / g(itail) * Egamma
	  endif

! ... CASE 3: phi computed in Friedrich prescription

	else if (extrad_flag.eq.3) then
	  if (itail.eq.0) stop 'Idiot! a multiplicative factor EXTRAD_PHI is not defined for peaking approx and EXTRAD_FLAG>2!'
	  if (itail.eq.1.or.itail.eq.2) then
	    x = Egamma/E(itail)
	    t = bt(itail)/etatzai
	    extrad_phi = extrad_phi * (1. - x + x**2/etatzai) *
     >		exp(t*((etatzai-0.5)-etatzai*x+x**2/2.)) * gamma(1.+bt(itail))
	  endif
	endif

	return
	end

!--------------------------------------------------------------------------
	real*8 function spen(x)

	implicit none

	integer	i
	real*8	x,s,y

! Approximation formula for the Spence-function or -integral according to
! abramowitz and stegun : formula 27.7.2

	y = 1.0
	s = 0.0
	i = 0
	do while (i.le.100 .and. abs(y).gt.abs(s)*1.d-4)
	  i = i+1
	  y = x*y
	  s = s+y/float(i**2)
	enddo
	spen = s

	return
	end

