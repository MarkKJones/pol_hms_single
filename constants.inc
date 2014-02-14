! CONSTANTS.INC

! Note: All angles are in radians
!	All distances are in cm
!	All energies (momenta,masses) are in MeV
!	All deltas are in percent
!	All densities (thicknesses) are in g/cm3 (g/cm2)
!	All B fields are in kG
!	All roads lead to infinity

! ... some physical constants
	real*8 me, Me2, Mp, Mp2, Mn, Mn2, Mpi, Mpi2, Mk, Mk2, Md, Md2
	real*8 Mlambda, Msigma0, Msigma_minus
	real*8 amu, hbarc, pi, twopi, alpha, alpi, degrad, euler
	real*8 Q_E, N_A
	parameter (me=0.51099906d00)
	parameter (Me2=Me**2)
	parameter (Mp=938.27231d00)
	parameter (Mp2=Mp**2)
	parameter (Mn=939.56563d00)
	parameter (Mn2=Mn**2)
	parameter (Mpi=139.56995d00)
	parameter (Mpi2=Mpi**2)
        parameter (Mk=493.677d00)
        parameter (Mk2=Mk**2)
        parameter (Md=1875.613d00)
        parameter (Md2=Md**2)
	parameter (Mlambda=1115.68d00)
	parameter (Msigma0=1192.64d00)
	parameter (Msigma_minus=1197.45d00)
	parameter (amu=931.49432d00)
	parameter (hbarc=197.327053d00)
	parameter (pi=3.141592653589793d00)
	parameter (twopi=2.*pi)
	parameter (alpha=1./137.0359895d00)
	parameter (alpi=alpha/pi)
	parameter (degrad=180.d0/pi)
	parameter (euler=0.577215665d00)
        PARAMETER (Q_E = 1.602d00)            !e- charge in uCoul (*1E-13)
        PARAMETER (N_A = 6.022d00)            !Avogadro's number (*1E+23)
