      SUBROUTINE QUASIY8_2006(Z,A,E0,EP,TH,SIGMA) 
      implicit none
      real*8 Z, A, E0, ep, th, sigma, F1, F2
      real*8 pi/3.1415927/, am/0.9383/, CHBAR/19732.8/
      real*8 ALPHA/7.29735E-03/
      real*8 sin2, cos2, tan2, q2, nu, wsq, W1, W2

      sin2 = sin(th * pi / 360.)**2
      cos2 = 1. - sin2
      tan2 = sin2 / cos2
      q2 = 4. * e0 * ep * sin2
      nu = e0 - ep
      wsq = am**2 + 2. * am * nu - q2
      call f1f2qe06(z, a, q2, wsq, f1, f2)
      W2 = F2 / nu
      W1 = F1 / am
      sigma = CHBAR**2 * 0.001 * ALPHA**2 / 4. * cos2 / 
     >   (e0 * sin2)**2 * (W2 + 2. * tan2 * W1)
      return
      end


      SUBROUTINE F1F2QE06(Z, A, qsq, wsq, F1, F2) 
c
C Calculates quasielastic A(e,e')X structure functions F1 and F2 PER NUCLEUS
c for A>2 uses superscaling from Sick, Donnelly, Maieron, nucl-th/0109032
c for A=2 uses pre-integrated Paris wave function (see ~bosted/smear.f)
c coded by P. Bosted August to October, 2006
c
c input: Z, A  (real*8) Z and A of nucleus (shoud be 2.0D0 for deueron)
c        Qsq (real*8) is 4-vector momentum transfer squared (positive in
c                     chosen metric)
c        Wsq (real*8) is invarinat mass squared of final state calculated
c                     assuming electron scattered from a free proton
c                 
c outputs: F1, F2 (real*8) are structure functions per nucleus
c
c Note: Deuteron agrees well with Laget (see ~bosted/eg1b/laget.f) for
c a) Q2<1 gev**2 and dsig > 1% of peak: doesnt describe tail at high W
c b) Q2>1 gev**2 on wings of q.e. peak. But, this model is up
c    to 50% too big at top of q.e. peak. BUT, F2 DOES agree very
c    nicely with Osipenko et al data from CLAS, up to 5 GeV**2

      IMPLICIT NONE     
      REAL*8 Z, A, avgN, F1, F2, wsq, qsq
      REAL*8 amp/0.93828/, amd/1.8756/
      REAL*8 PAULI_SUP1, PAULI_SUP2
      REAL*8 GEP, GEN, GMP, GMN, Q, Q3, Q4
      REAL*8 RMUP/ 2.792782/ ,RMUN/ -1.913148 /                
      real*8 pz, nu, dpz, pznom, pzmin
      real*8 qv, TAU, W1, W2, FY, dwmin, w2p
      real kappa, lam, lamp, taup, squigglef, psi, psip, nuL, nut
      real kf, es, GM2bar, GE2bar, W1bar, W2bar, Delta, GL, GT
      integer IA, izz, izzmin, izp, izznom, izdif

c Look up tables for deuteron case
       real*8 fyd(200)/
     > 0.00001,0.00002,0.00003,0.00005,0.00006,0.00009,0.00010,0.00013,
     > 0.00015,0.00019,0.00021,0.00026,0.00029,0.00034,0.00038,0.00044,
     > 0.00049,0.00057,0.00062,0.00071,0.00078,0.00089,0.00097,0.00109,
     > 0.00119,0.00134,0.00146,0.00161,0.00176,0.00195,0.00211,0.00232,
     > 0.00252,0.00276,0.00299,0.00326,0.00352,0.00383,0.00412,0.00447,
     > 0.00482,0.00521,0.00560,0.00603,0.00648,0.00698,0.00747,0.00803,
     > 0.00859,0.00921,0.00985,0.01056,0.01126,0.01205,0.01286,0.01376,
     > 0.01467,0.01569,0.01671,0.01793,0.01912,0.02049,0.02196,0.02356,
     > 0.02525,0.02723,0.02939,0.03179,0.03453,0.03764,0.04116,0.04533,
     > 0.05004,0.05565,0.06232,0.07015,0.07965,0.09093,0.10486,0.12185,
     > 0.14268,0.16860,0.20074,0.24129,0.29201,0.35713,0.44012,0.54757,
     > 0.68665,0.86965,1.11199,1.43242,1.86532,2.44703,3.22681,4.24972,
     > 5.54382,7.04016,8.48123,9.40627,9.40627,8.48123,7.04016,5.54382,
     > 4.24972,3.22681,2.44703,1.86532,1.43242,1.11199,0.86965,0.68665,
     > 0.54757,0.44012,0.35713,0.29201,0.24129,0.20074,0.16860,0.14268,
     > 0.12185,0.10486,0.09093,0.07965,0.07015,0.06232,0.05565,0.05004,
     > 0.04533,0.04116,0.03764,0.03453,0.03179,0.02939,0.02723,0.02525,
     > 0.02356,0.02196,0.02049,0.01912,0.01793,0.01671,0.01569,0.01467,
     > 0.01376,0.01286,0.01205,0.01126,0.01056,0.00985,0.00921,0.00859,
     > 0.00803,0.00747,0.00698,0.00648,0.00603,0.00560,0.00521,0.00482,
     > 0.00447,0.00412,0.00383,0.00352,0.00326,0.00299,0.00276,0.00252,
     > 0.00232,0.00211,0.00195,0.00176,0.00161,0.00146,0.00134,0.00119,
     > 0.00109,0.00097,0.00089,0.00078,0.00071,0.00062,0.00057,0.00049,
     > 0.00044,0.00038,0.00034,0.00029,0.00026,0.00021,0.00019,0.00015,
     > 0.00013,0.00010,0.00009,0.00006,0.00005,0.00003,0.00002,0.00001/
       real*8 avp2(200)/
     >     1.0,0.98974,0.96975,0.96768,0.94782,0.94450,0.92494,0.92047,
     > 0.90090,0.89563,0.87644,0.87018,0.85145,0.84434,0.82593,0.81841,
     > 0.80021,0.79212,0.77444,0.76553,0.74866,0.73945,0.72264,0.71343,
     > 0.69703,0.68740,0.67149,0.66182,0.64631,0.63630,0.62125,0.61154,
     > 0.59671,0.58686,0.57241,0.56283,0.54866,0.53889,0.52528,0.51581,
     > 0.50236,0.49291,0.47997,0.47063,0.45803,0.44867,0.43665,0.42744,
     > 0.41554,0.40656,0.39511,0.38589,0.37488,0.36611,0.35516,0.34647,
     > 0.33571,0.32704,0.31656,0.30783,0.29741,0.28870,0.27820,0.26945,
     > 0.25898,0.25010,0.23945,0.23023,0.21943,0.20999,0.19891,0.18911,
     > 0.17795,0.16793,0.15669,0.14667,0.13553,0.12569,0.11504,0.10550,
     > 0.09557,0.08674,0.07774,0.06974,0.06184,0.05484,0.04802,0.04203,
     > 0.03629,0.03129,0.02654,0.02247,0.01867,0.01545,0.01251,0.01015,
     > 0.00810,0.00664,0.00541,0.00512,0.00512,0.00541,0.00664,0.00810,
     > 0.01015,0.01251,0.01545,0.01867,0.02247,0.02654,0.03129,0.03629,
     > 0.04203,0.04802,0.05484,0.06184,0.06974,0.07774,0.08674,0.09557,
     > 0.10550,0.11504,0.12569,0.13553,0.14667,0.15669,0.16793,0.17795,
     > 0.18911,0.19891,0.20999,0.21943,0.23023,0.23945,0.25010,0.25898,
     > 0.26945,0.27820,0.28870,0.29741,0.30783,0.31656,0.32704,0.33571,
     > 0.34647,0.35516,0.36611,0.37488,0.38589,0.39511,0.40656,0.41554,
     > 0.42744,0.43665,0.44867,0.45803,0.47063,0.47997,0.49291,0.50236,
     > 0.51581,0.52528,0.53889,0.54866,0.56283,0.57241,0.58686,0.59671,
     > 0.61154,0.62125,0.63630,0.64631,0.66182,0.67149,0.68740,0.69703,
     > 0.71343,0.72264,0.73945,0.74866,0.76553,0.77444,0.79212,0.80021,
     > 0.81841,0.82593,0.84434,0.85145,0.87018,0.87644,0.89563,0.90090,
     > 0.92047,0.92494,0.94450,0.94782,0.96768,0.96975,0.98974,1.0/

! return if proton: future change this to allow for
! equivalent W resolution
      F1 = 0.
      F2 = 0.
      IA = int(A)
      avgN = A - Z 
      IF (iA.EQ.1) RETURN                                      

! some kinematic factors. Return if nu or qsq is negative
      Nu = (wsq - amp**2 + qsq) / 2. / amp
      if(nu .le. 0.0 .or. qsq .lt. 0.) return
      TAU   = QSQ / 4.0 / amp**2                                        
      qv = sqrt(nu**2 + qsq)

! Bosted fit for nucleon form factors Phys. Rev. C 51, p. 409 (1995)
      Q = sqrt(QSQ)
      Q3 = QSQ * Q
      Q4 = QSQ**2
      GEP = 1./  (1. + 0.14 * Q + 3.01 * QSQ + 0.02 * Q3 + 
     >  1.20 * Q4 + 0.32 * Q**5)
      GMP = RMUP * GEP
      GMN = RMUN / (1.- 1.74 * Q + 9.29 * QSQ - 7.63 * Q3 + 
     >  4.63 * Q4)
      GEN = 1.25 * RMUN * TAU / (1. + 18.3 * TAU) / 
     >  (1. + QSQ / 0.71)**2

! Get kf and Es from superscaling from Sick, Donnelly, Maieron,
c nucl-th/0109032
      if(IA.eq.2) kf=0.085
      if(iA.eq.2) Es=0.0022
      if(IA.eq.3) kf=0.150
      if(iA.eq.3) Es=0.010 
      if(IA.eq.4) kf=0.200
      if(iA.eq.4) Es=0.015 
      if(IA.gt.4) kf=0.165
      if(iA.gt.4) Es=0.015 
      if(IA.gt.7) kf=0.228
      if(iA.gt.7) Es=0.020 
      if(IA.gt.16) kf=0.230
      if(iA.gt.16) Es=0.025 
      if(IA.gt.25) kf=0.236
      if(iA.gt.25) Es=0.018 
      if(IA.gt.38) kf=0.241
      if(iA.gt.38) Es=0.028 
      if(IA.gt.55) kf=0.241
      if(iA.gt.55) Es=0.023 
      if(IA.gt.60) kf=0.245
      if(iA.gt.60) Es=0.028 

! Pauli suppression model from Tsai RMP 46,816(74) eq.B54
      IF((QV .GT. 2.* kf).OR.(iA.EQ.1)) THEN
        PAULI_SUP2 =1.0
      ELSE
        PAULI_SUP2 = 0.75 * (QV / kf) * (1.0 - ((QV / kf)**2)/12.)
      ENDIF
      PAULI_SUP1 = PAULI_SUP2

! structure functions with off shell factors
      kappa = qv / 2. / amp
      lam = nu / 2. / amp
      lamp = lam - Es / 2. / amp
      taup = kappa**2 - lamp**2
      squigglef = sqrt(1. + (kf/amp)**2) -1.
! Very close to treshold, could have a problem
      if(1.+lamp.le.0.) return
      if(taup * (1. + taup).le.0.) return

      psi =  (lam  - tau ) / sqrt(squigglef) /
     >  sqrt((1.+lam )* tau + kappa * sqrt(tau * (1. + tau)))

      psip = (lamp - taup) / sqrt(squigglef) / 
     >  sqrt((1.+lamp)*taup + kappa * sqrt(taup * (1. + taup)))

      nuL = (tau / kappa**2)**2

c changed definition of nuT from
c      nuT = tau / 2. / kappa**2 + tan(thr/2.)**2
c to this, in order to separate out F1 and F2 (F1 prop. to tan2 term)
      nuT = tau / 2. / kappa**2 

      GM2bar = Pauli_sup1 * (Z * GMP**2 + avgN * GMN**2)  
      GE2bar = Pauli_sup2 * (Z * GEP**2 + avgN * GEN**2) 
      W1bar = tau * GM2bar
      W2bar = (GE2bar + tau * GM2bar) / (1. + tau)

      Delta = squigglef * (1. - psi**2) * (
     >  sqrt(tau * (1.+tau)) / kappa + squigglef/3. *
     >  (1. - psi**2) * tau / kappa**2)

      GL = kappa**2 / tau * (GE2bar + Delta * W2bar) / 
     >  2. / kappa / (1. + squigglef * (1. + psi**2) / 2.)
      GT = (2. * tau * GM2bar + Delta * W2bar) /
     >  2. / kappa / (1. + squigglef * (1. + psi**2) / 2.)

! from Maria Barbaro: see Amaro et al., PRC71,015501(2005).
      FY = 1.5576 / (1. + 1.7720**2 * (psip + 0.3014)**2) / 
     >   (1. + exp(-2.4291 * psip)) / kf

! Use PWIA and Paris W.F. for deuteron to get better FY
      if(IA.eq.2) then
! value assuming average p2=0.
        pz = (qsq - 2. * amp * nu ) / 2. / qv
        izz = int((pz + 1.0) / 0.01) + 1
        izz = min(200,max(1,izz))
        izznom = izz
! ignoring energy term, estimate change in pz to compensate
! for avp2 term
        dpz = avp2(izznom) / 2. / qv
        izdif = dpz * 150. 
        dwmin=1.E6
        izzmin=0
        do izp = izznom, min(200, max(1, izznom + izdif))
          pz = -1. + 0.01 * (izp-0.5)
c *** this version gives worse agreement with laget than
c         w2p = (amd + nu - sqrt(am**2 + avp2(izp)))**2 - 
c    >      qv**2 + 2. * qv * pz - avp2(izp)
c this version!
          w2p = (amd + nu - amp )**2 - 
     >      qv**2 + 2. * qv * pz - avp2(izp)
c if passed first minimum, quit looking so don't find second one
          if(abs(w2p - amp**2).gt.dwmin) goto 11
          if(abs(w2p - amp**2).lt.dwmin) then
            dwmin = abs(w2p - amp**2)
            izzmin = izp
          endif
        enddo
 11     izz = min(199,max(2,izzmin))
! search for minimum in 1/10th bins locally
        pznom = -1. + 0.01 * (izz-0.5)
        dwmin=1.E6
        do izp = 1,19
          pz = pznom - 0.01 + 0.001 * izp
          w2p = (amd + nu - amp )**2 - 
     >      qv**2 + 2. * qv * pz - avp2(izz)
          if(abs(w2p - amp**2).lt.dwmin) then
            dwmin = abs(w2p - amp**2)
            pzmin = pz
          endif
        enddo
        if(dwmin.ge.1.e6.or.abs(pznom-pzmin).gt.0.01) 
     >     write(6,'(1x,''error in dwmin,pzmin'')')
        if(pzmin.lt.pznom) then
          fy = fyd(izz) - (fyd(izz-1) - fyd(izz)) * 
     >      (pzmin - pznom) / 0.01
        else
          fy = fyd(izz) + (fyd(izz+1) - fyd(izz)) * 
     >      (pzmin - pznom) / 0.01
        endif
      endif

c final results
      F2 = nu * FY * (nuL * GL + nuT * GT)
      F1 = amp * FY * GT / 2.

      return
      end
