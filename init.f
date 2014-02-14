!------------------------------------------------------------------------
	subroutine radc_init(Z)

	implicit none
c	include 'simulate.inc'
	include 'radc.inc'
	include 'brem.inc'
        real*8 L1, L2, Z
!------------------------------------------------------------------------
! First, about those (mysterious) 2 main 'radiative option' flags ...
!
! The significance of RAD_FLAG:
!	 RAD_FLAG  = 0	.. use best available formulas, generate in
!			.. (ntail,Egamma) basis
!		   = 1	.. use BASICRAD only, generate in (ntail,Egamma)
!			.. basis
!		   = 2	.. use BASICRAD only, generate in (Egamma(1,2,3))
!			.. basis but prevent overlap of tails (bogus, note)
!		   = 3	.. use BASICRAD only, generate in (Egamma(1,2,3))
!			.. allowing radiation in all 3 directions
!			.. simultaneously
! The (ntail,Egamma) basis can be called the PEAKED basis since it allows
! only 3 photon directions. PEAKED_BASIS_FLAG is set to zero below when
! the peaked basis is being used, in this way we can conveniently tell
! the BASICRAD routine to use the full Egamma formula to generate the gamma
! energy whenever it's called.
!
! (See N. Makins' thesis, section 4.5.6, for more help on this point)
!
! The significance of EXTRAD_FLAG:
!	EXTRAD_FLAG = 1	.. use only BASIC external radiation formulas
!			.. (phi = 1)
!		    = 2	.. use BASIC ext rad formulas x phi
!		    = 3 .. use Friedrich approximation the way we always
!			.. have
!		    = 0 .. use DEFAULTS: 3 for RAD_FLAG = 0, 1 otherwise; note
!			   that the defaults mimic the hardwired 'settings'
!			   in SIMULATE, which doesnt read EXTRAD_FLAG
!			   but determines what to do based on RAD_FLAG
!-------------------------------------------------------------------------

! Check setting of EXTRAD_FLAG

        rad_flag = 0
        extrad_flag = 2
c	if (extrad_flag.eq.0) then
c	  if (rad_flag.eq.0) then
c	    extrad_flag = 3
c	  else if (rad_flag.eq.1 .or. rad_flag.eq.2 .or. rad_flag.eq.3) then
c	    extrad_flag = 1
c	  endif
c	else if (extrad_flag.lt.0) then
c	  stop 'Imbecile! check your stupid setting of EXTRAD_FLAG'
c	endif

! 'etatzai' parameter

	L1 = log(184.15) - log(Z)/3.0
	L2 = log(1194.) - 2.*log(Z)/3.0
	if(Z.eq.1)then
	  L1=5.31
	  L2=6.144
	endif
	etatzai = (12.0+(Z+1.)/(Z*L1+L2))/9.0

! Initialize brem flags (brem doesn't include the normal common blocks)

c	exponentiate = use_expon
	exponentiate = .false.
	include_hard = .true.
	calculate_spence = .true.

	return
	end

!---------------------------------------------------------------------

	subroutine radc_init_ev (tb,ta,Ein,Ep,th_ev,ue)

	implicit none
c	include 'structures.inc'
	include 'radc.inc'

	integer	i
        real*8 tb,ta,Ein,Ep,th_ev,ue(3)
	real*8 r, Ecutoff, dsoft, dhard, dsoft_prime
	real*8 lambda_dave, bremos
	real*8 zero
	parameter (zero=0.0d0)	!double precision zero for subroutine calls.

! Compute some quantities that will be needed for rad corr on this event

! ... factor for limiting energy of external radiation along incident electron
!	etta = 1.0 + 2*vertex.ein*sin(vertex.e.theta/2.)**2/(targ.A*amu)
! ... moron move! let's can that etta factor ...

	etta = 1.0

! ... the bt's

	bt(1) = etatzai*tb     !length before vertex in rad length
        bt(2) = etatzai*ta     !length after vertex in rad length

! ... the lambda's (effective bt's for internal radiation)

c	doing_tail(1) = one_tail.eq.0 .or. one_tail.eq.1 
c	doing_tail(2) = one_tail.eq.0 .or. one_tail.eq.2

        doing_tail(1) = .true.
        doing_tail(2) = .true.
        doing_tail(3) = .false.

	do i=1,2
	  lambda(i) = lambda_dave(i,1,doing_tail(3),Ein,Ep,th_ev)
	enddo

! ... get the hard correction factor. don't care about Ecutoff! Just want dhard here

	Ecutoff = 450.
	r = bremos(Ecutoff, zero, zero, Ein, Ep*ue(1), Ep*ue(2), Ep*ue(3), 
     >             dsoft, dhard, dsoft_prime)
	hardcorfac = 1./(1.-dhard)
	g(4)=-dsoft_prime*Ecutoff+bt(1)+bt(2)

! ... initialize the parameters needed for our "basic" calculation

	call basicrad_init_ev (Ein,Ep)

! ... the relative magnitudes of the three tails (we may not need them)

	do i=1,2
	  frac(i) = g(i)/g(0)
	enddo

	return
	end

!--------------------------------------------------------------------------
	real*8 function lambda_dave(itail,plus_flag,doing_proton,e1,e2,th)

	implicit none
	include 'constants.inc'

	integer		itail,plus_flag
	real*8		e1,e2,th
	real*8		plus_term
	logical		doing_proton, warned/.false./

! The extended peaking approximation

	plus_term = 0.0
	if (plus_flag.eq.1 .and. itail.lt.3) then
	  plus_term = log((1.-cos(th))/2.)

! ... only add in term due to ep interference if we're using proton
! ... radiation

	  if (doing_proton) plus_term = plus_term + 2.*log(e1/e2)
	endif

! Compute lambdas

	if (itail.eq.1) then
	  lambda_dave = alpi*(2.*log(2.*e1/Me) -1. + plus_term)
	else if (itail.eq.2) then  
	  lambda_dave = alpi*(2.*log(2.*e2/Me) -1. + plus_term)
	endif

	return
	end

!-----------------------------------------------------------------------
	subroutine basicrad_init_ev (e1,e2)

	implicit none
c	include 'simulate.inc'
	include 'radc.inc'
	include 'constants.inc'

	real*8 one
	parameter (one=1.)

	integer i
	real*8 e1,e2,e(2),gamma

	e(1) = e1
	e(2) = e2

! bt's for internal + external
! ??? One possibility for shutting off tails 1 or
! 2 is to set g(1/2) = 0 here ... note that lambda(3) is set to 0 in
! lambda_dave at the moment, AND proton terms are removed from brem if proton
! radiation off. something analogous and similarly consistent would have to be
! done for the other tails, right now they're just nixed in generate_rad. also
! check ALL ntail.eq.0 checks in kinema constraint lines of generate_rad

	g(1) = lambda(1) + bt(1)
	g(2) = lambda(2) + bt(2)
	g(0) = g(1)+g(2)

! Internal constants

	c_int(1) = lambda(1)/(e(1)*e(2))**(lambda(1)/2.)
	c_int(2) = lambda(2)/(e(1)*e(2))**(lambda(2)/2.)

	do i = 1, 2
	  c_int(i) = c_int(i) * exp(-euler*lambda(i)) / gamma(one+lambda(i))
	enddo

	g_int = lambda(1) + lambda(2) 
	c_int(0) = c_int(1)*c_int(2) * g_int / lambda(1)/lambda(2)
	c_int(0) = c_int(0) * gamma(one+lambda(1)) * gamma(one+lambda(2))
     >             / gamma(one+g_int)

! External constants

	do i = 1, 2
	  c_ext(i) = bt(i)/e(i)**bt(i)/gamma(one+bt(i))
	enddo
	g_ext = bt(1) + bt(2)
	c_ext(0) = c_ext(1)*c_ext(2) * g_ext / bt(1)/bt(2)
	c_ext(0) = c_ext(0)*gamma(one+bt(1))*gamma(one+bt(2))/gamma(one+g_ext)

! Internal + external constants

	do i = 1, 2
	  c(i) = c_int(i) * c_ext(i) * g(i)/lambda(i)/bt(i)
     >		* gamma(one+lambda(i))*gamma(one+bt(i))/gamma(one+g(i))
	enddo

! Finally, constant for combined tails

	c(0) = c(1)*c(2) * g(0)/g(1)/g(2)
	c(0)=c(0)*gamma(one+g(1))*gamma(one+g(2))/gamma(one+g(0))
	c(4)=g(4)/(e1*e2)**g(4)/gamma(one+g(4))

	return
	end

!-------------------------------------------------------------------------
	real*8 function gamma(x)

	implicit none

	integer i, n, s
	real*8	x, y

! Compute gamma function for xin
! The series computes gamma(1+y), and must be fed y between 0 and 1
! This can be accomplished using the relation
! gamma(y+n+1) = (y+n)*gamma(y+n) = ... = (y+n)*...*(y+1)*gamma(y+1)

	gamma = 1.0
	n = nint((x-1)-0.5)
	y = x-1 - n
	if (n.ne.0) then
	  s = sign(1,n)
	  do i = s, n, s
	    gamma = gamma*(y+1+i)**s
	  enddo
	endif
	gamma = gamma * (1. - 0.5748646*y + 0.9512363*y**2 -
     >		 0.6998588*y**3 + 0.4245549*y**4 - 0.1010678*y**5)
	return
	end


