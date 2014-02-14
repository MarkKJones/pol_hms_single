      subroutine pb_ext_sub(ZZ,AA,ebeam,nu,theta,sigma)
c
c       ebeam ( Mev), nu ( MeV), th (deg)
c       return sigma (nb/GeV/sr)
      Implicit None 
      COMMON /TARGT/ iZ,iA,avgN,avgA,avgM,amuM 
      integer iz,ia
      REAL avgN, avgA, avgM, amuM
      real*8 ebeam,nu,theta              
      real e0,ep,th,sigma_inel,sigma_qf
      real*8 sigma,ZZ,AA
      COMMON/IKK12/IG,IDUT,INEL_MODEL,PAULI_MODEL,NUC_METHOD,NUC_MODEL
      INTEGER IG,IDUT,INEL_MODEL,PAULI_MODEL,NUC_METHOD,NUC_MODEL
c 
      ia=int(aa)
      iz=int(zz)
      ig=15
      idut=13
      inel_model = 0
      pauli_model = 1  
      nuc_method = 1
      nuc_model =1
c
      call weiz ! calculates avgM, sets avgA=iA
      e0=ebeam/1000.
      ep=(ebeam-nu)/1000.
      th=theta
      avgN = avgA-iZ                                                    
      amuM = avgM/.931501                                               
      call secnuclw(e0,ep,th,sigma_inel)
      CALL  QUASIY8(E0,EP,TH,SIGMA_qf)
      sigma = sigma_inel + sigma_qf
      return
c
 999  end

                                                                       
