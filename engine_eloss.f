********************************************************************************

      subroutine loss(particle,z,a,thick,dens,velocity,e_loss)
*-------------------------------------------------------------
*- 	 Prototype C function 
*- 
*-
*-  Purpose and Method :  Calculate energy loss 
*-  
*-  Output: -
*-  Created   1-Dec-1995  Rolf Ent
*- 
*-  Verification:  The non-electron portion on this subr. is Bethe_Bloch
*- 		   equation (Physial Review D vol.50 (1994) 1251 with full
*-      	   calculation of Tmax and the density correction. The electron
*-      	   part has been switched from O'Brien, Phys. Rev. C9(1974)1418,
*-      	   to Bethe-Bloch with relativistic corrections and density
*-      	   density correction, Leo, Techniques for Nuclear and Particle 
*-      	   Physics Experiments
*- 		   J. Volmer 8/2/98 16:50
*------------------------------------------------------------------------------*

      IMPLICIT NONE
      SAVE

      INTEGER*4 particle        ! 0 e; 1 p; 2 neutron
      REAL*8 eloss,z,a
      REAL*8 thick       ! grams/cm**2
      REAL*8 dens        ! grams/cm**3
      REAL*8 beta,e_loss
      REAL*8 icon_ev
      REAL*8 icon_gev
      REAL*8 denscorr,hnup,c0,log10bg,pmass,tmax,gamma,velocity
      REAL*8 tau,betagamma

      real*8 npartmass
      real*8 mass_electron
      real*8 hpartmass
      real*8 mass_nucleon
      integer*4 ELECTRON,PROTON,NEUTRON
      parameter (ELECTRON = 0)
      parameter (PROTON   = 1)
      parameter (NEUTRON  = 2)

      LOGICAL ElossDebug
      parameter (ElossDebug = .false.)
 
      mass_electron = .511d-3
      hpartmass = .93827d-3
      mass_nucleon =  .93827d-3
      e_loss = 0.0
      eloss  = 0.0
*     csa -- This should be an externally supplied value, but
*     it is not yet in the NDET code
      if (particle.eq.PROTON) then
         npartmass      = mass_nucleon  !proton=0.93827
      else
         npartmass      = 0.93955 !neutron
      endif


*****************************************************************************
* calculate the mean excitation potential I in a newer parametrization 
* given in W.R. Leo's Techniques for Nuclear and Particle Physics Experiments
*****************************************************************************

      if (z.lt.1.5) then
         icon_ev = 21.8
      elseif (z.lt.13) then
         icon_ev = 12.*z+7.
      else
         icon_ev = z*(9.76+58.8*z**(-1.19))
      endif
      icon_gev = icon_ev*1.0e-9 ! from 9
**********************************************
* extract the velocity of the particle:
*     hadrons:   velocity = beta
*     electrons: velocity = log_10(beta*gamma)
**********************************************

      if (particle.eq.ELECTRON) then
         log10bg=velocity
         betagamma=exp(velocity*log(10.))
         beta=betagamma/(sqrt(1.+betagamma**2))
         gamma=sqrt(1.+betagamma**2)
         tau=gamma-1.

      elseif (particle.eq.PROTON) then
         beta=abs(velocity)
         if (beta.ge.1.) beta=.9995
         if (beta.le.0.1) beta=0.1
         gamma=1./sqrt(1.-beta**2)
         betagamma=beta*gamma
         log10bg=log(betagamma)/log(10.)
         tau=gamma-1.

      elseif (particle.eq.NEUTRON) then     ! ignored for now!
         e_loss = 0.

      else
         write(6,*)' HEY -- bogus particle type!'
         RETURN
      endif

******************************************************
* calculate the density correction, as given in Leo,
* with Sternheimer's parametrization
* I is the mean excitation potential of the material
* hnup= h*nu_p is the plasma frequency of the material
******************************************************

      if(A.gt.0.) then
         HNUP=28.816E-9*sqrt(abs(DENS*Z/A)) ! from 9
      else
         HNUP=28.816E-9*sqrt(abs(DENS*Z/1.))
      endif

      C0 = -2.*(log(icon_gev)-log(hnup)+.5)   ! log(a/b) = log(a) - log(b)

      if(log10bg.lt.0.) then
         denscorr=0.
      elseif(log10bg.lt.3.) then
         denscorr=C0+2*log(10.)*log10bg+abs(C0/27.)*(3.-log10bg)**3
      elseif(log10bg.lt.4.7) then
         denscorr=C0+2*log(10.)*log10bg
      else
         denscorr=C0+2*log(10.)*4.7
      endif


**********************************************************************       
* now calculate the energy loss for electrons 
**********************************************************************

      if (particle.eq.ELECTRON) then
         if((thick.gt.0.0).and.(dens.gt.0.0).and.(a.gt.0.).and.(beta.gt.0.)
     >       .and.(tau.gt.0).and.(betagamma.gt.0))then
            eloss=0.1535e-03*z/a*thick/beta**2*(
     >           2*log(tau)+log((tau+2.)/2.)
     >            -2*(log(icon_gev)-log(mass_electron))
     >           +1-beta**2+(tau**2/8-(2*tau+1)*log(2.))/(tau+1)**2
     >           -(-(2*(log(icon_gev)-log(hnup))+1)+2*log(betagamma)))
         endif

********************************************************************      
* now calculate the energy loss for hadrons 
********************************************************************

      elseif(particle.eq.PROTON) then

*        first calculate the maximum possible energy transfer to
*        an orbital electron, find out what the hadron mass is
         pmass=max(hpartmass,npartmass)
         if (pmass.lt.2*mass_electron) pmass=0.5
         tmax=abs(2*mass_electron*beta**2*gamma**2/
     >        (1+2*abs(gamma)*mass_electron/pmass+(mass_electron/pmass)**2))

         if((thick.gt.0.0).and.(beta.gt.0.0)
     >            .and.(beta.lt.1.0).and.(a.gt.0.))then

            eloss = abs(2.*0.1535e-3*Z/A*thick/beta**2)*
     >           ( .5*(log(2*mass_electron) + 2*log(beta) + 2*log(gamma)
     >                 + log(tmax) - 2*log(icon_gev))
     >            -beta**2-denscorr/2.)

         endif
	 
      elseif (particle.eq.NEUTRON) then   ! ignored for now!
         e_loss = 0.
	 
      else
         e_loss = 0.
      endif

      if (eloss.le.0.) write(6,*)'loss: eloss<=0!'

      e_loss = eloss*1000    ! units should be in MeV


      if ((ElossDebug).or.(eloss.le.0)) then
         write(6,*) '************************************************'
         write(6,91) 'particle','ztgt','atgt','thick','dens','velocity','e_loss'
         write(6,90) particle,z,a,thick,dens,velocity,e_loss
         write(6,*) ' '
         write(6,'(4A10)') 'velocity','beta','pmass','denscorr'
         write(6,'(6(2x,f8.5))') velocity,beta,pmass,denscorr
         write(6,'(6A10)') 'betagamma','log10bg',
     >                                'tau','gamma','icon_ev','hnup (eV)'
         write(6,'(6(2x,F8.3))') betagamma,log10bg,tau,gamma,icon_ev,hnup*1e9
         write(6,*) ' '
      endif

 91   format('loss: ',7(A11))
 90   format('loss: ',2x,I9,7(2x,f10.6))


      RETURN
      END
