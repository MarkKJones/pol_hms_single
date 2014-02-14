CCC  Version 072205  -  Author:  M.E. Christy                        CCC
CCC  Fit form is empirical.  Please do not try to interpret physics  CCC
CCC  from it.  This routine needs the parameter files f1parms.dat    CCC
CCC  and flparms.dat.  Units are ub/Sr/Gev.                          CCC

             
      SUBROUTINE RESMOD(sf,w2,q2,xval,sig) 

      IMPLICIT NONE
      REAL*8 W,w2,q2,mp,mp2,xb,xth(4),sig,xval(50),mass(7),width(7)
      REAL*8 height(7),sig_del,sig_21,sig_22,sig_31,sig_32,rescoef(6,4)
      REAL*8 nr_coef(3,4),wdif,wdif2,wr,sig_nr,sig_4,q2temp
      REAL*8 mpi,meta,intwidth(6),k,kcm,kcmr(6),ppicm,ppi2cm,petacm
      REAL*8 ppicmr(6),ppi2cmr(6),petacmr(6),epicmr(6),epi2cmr(6)
      REAL*8 eetacmr(6),epicm,epi2cm,eetacm,br_21_1,br_21_2
      REAL*8 sig_res,sig_4L,sigtemp,slope,q2low
      INTEGER i,j,l,lmax,num,sf
      LOGICAL lowq2

      lowq2 = .false.
      lmax = 1
      q2temp = q2

      mp = 0.9382727
      mpi = 0.136
      meta = 0.547
      mp2 = mp*mp
      W = sqrt(w2)
      wdif = w - (mp + mpi)
      wr = wdif/w


      if(sf.EQ.1) then
        q2low = 0.15
      else
        q2low = 0.05
      endif 

      if(q2.LT.q2low) then
        lowq2 = .true.
        lmax = 2
      endif

c      write(6,*) q2,lmax

      do l=1,lmax

               
        if(l.EQ.1.AND.lowq2) then
          q2 = q2low
        elseif(l.EQ.2.AND.lowq2) then
          q2 = q2low + 0.1
        endif

        xb = q2/(q2+w2-mp2)
        xth(1) = (q2 + xval(50))/(w2-mp2-0.136+q2)


CCC  Calculate kinematics needed for threshold Relativistic B-W  CCC

        k = (w2 - mp2)/2./mp
        kcm = (w2-mp2)/2./w

        epicm = (W2 + mpi**2 -mp2 )/2./w
        ppicm = SQRT(MAX(0.0,(epicm**2 - mpi**2)))
        epi2cm = (W2 + (2.*mpi)**2 -mp2 )/2./w
        ppi2cm = SQRT(MAX(0.0,(epi2cm**2 - (2.*mpi)**2)))
        eetacm = (W2 + meta*meta -mp2 )/2./w
        petacm =  SQRT(MAX(0.0,(eetacm**2 - meta**2)))

        num = 0

        do i=1,6
          num = num + 1
          mass(i) = xval(i)

          kcmr(i) = (mass(i)**2.-mp2)/2./mass(i)
          epicmr(i) = (mass(i)**2 + mpi**2 -mp2 )/2./mass(i)
          ppicmr(i) = SQRT(MAX(0.0,(epicmr(i)**2 - mpi**2)))
          epi2cmr(i) = (mass(i)**2 + (2.*mpi)**2 -mp2 )/2./mass(i)
          ppi2cmr(i) = SQRT(MAX(0.0,(epi2cmr(i)**2 - (2.*mpi)**2)))
          eetacmr(i) = (mass(i)**2 + meta*meta -mp2 )/2./mass(i)
          petacmr(i) =  SQRT(MAX(0.0,(eetacmr(i)**2 - meta**2)))

        enddo
 
        do i=1,6
          num = num + 1
          intwidth(i) = xval(num)
          width(i) = intwidth(i)
        enddo

           
c      write(6,*) "1:  ",num

        do i=1,6
          do j=1,4
            num = num + 1
            rescoef(i,j)=xval(num)
          enddo
          if(sf.EQ.1) then

            height(i) = rescoef(i,1)/
     &        (1.+q2/rescoef(i,2))**(rescoef(i,3)+rescoef(i,4)*q2)

            if(i.EQ.1) height(i) = 3.0*
     &         (1.+rescoef(i,1)*q2/(1.+rescoef(i,2)*q2)/
     &         (1.+q2/rescoef(i,3))**rescoef(i,4))/(1.+q2/0.71)**2. 
            if(i.EQ.2) height(i) = 1.4*
     &         (1.+rescoef(i,1)*q2/(1.+rescoef(i,2)*q2)/
     &         (1.+q2/rescoef(i,3))**rescoef(i,4))/(1.+q2/0.71)**2. 
            if(i.EQ.5) height(i) = 0.3*
     &         (1.+rescoef(i,1)*q2/(1.+rescoef(i,2)*q2)/
     &         (1.+q2/rescoef(i,3))**rescoef(i,4))/(1.+q2/0.71)**2. 

          else
c            height(i) = rescoef(i,1)*
c     &            (1.+q2/rescoef(i,2))**(-1.*rescoef(i,3))
            height(i) = rescoef(i,1)*q2/(1.+rescoef(i,4)*q2)/
     &            (1.+q2/rescoef(i,2))**rescoef(i,3)     !!!  Galster Shape  !!!             
          endif
          if(height(i).LT.0) height(i) = 0. 

        enddo
     

        do i=1,3
          do j=1,4
            num = num + 1
            nr_coef(i,j)=xval(num)
          enddo
        enddo

        if(sf.EQ.2) then      !!!  Put in Roper  !!!
          mass(7) = xval(41)
          width(7) = xval(42)
          height(7) = xval(43)*q2/(1.+xval(44)*q2)/
     &            (1.+q2/xval(45))**xval(46)     !!!  Galster Shape  !!!
          sig_4L   = height(7)*width(7)/((W-mass(7))**2. 
     &               + 0.25*width(7)*width(7))  
        else
          mass(7) = xval(47)
          width(7) = xval(48)
          height(7) = xval(49)/(1.+q2/0.61)**3.    
          sig_4L   = height(7)*width(7)/((W-mass(7))**2. 
     &               + 0.25*width(7)*width(7))  
        endif


CC   Calculate Breit-Wigners for the 3 resonance regions   CC


        sig_32 = width(5)/((W-mass(5))**2. 
     &               + 0.25*width(5)*width(5))
        sig_4   = width(6)/((W-mass(6))**2. 
     &               + 0.25*width(6)*width(6))

        if(sf.EQ.1) then 
         br_21_1 = 0.5
         br_21_2 = 0.5
        else
          br_21_1 = 0.985
          br_21_2 = 1.-br_21_1
        endif

        width(1)=intwidth(1)*ppicm/ppicmr(1)
        width(2)=intwidth(2)*(br_21_1*ppicm/ppicmr(2)
     &            +br_21_2*petacm/petacmr(2))
        width(3)=intwidth(3)*(0.5*ppicm/ppicmr(3)+0.5*ppi2cm/ppi2cmr(3))
        width(4)=intwidth(4)*
     &                     (0.65*ppicm/ppicmr(4)+0.35*ppi2cm/ppi2cmr(4))

c      write(6,*) ppicm,ppicmr(3),petacm,petacmr(3),intwidth(3)

        sig_del = ppicm/kcm/((W2 - mass(1)**2.)**2. 
     &              + (mass(1)*width(1))**2.)
        sig_21 =  (0.5*ppicm+0.5*petacm)/kcm/
     &           ((W2 - mass(2)**2.)**2. + (mass(2)*width(2))**2.)
        sig_22 =  (0.5*ppicm+0.5*ppi2cm)/2./kcm/
     &           ((W2 - mass(3)**2.)**2. + (mass(3)*width(3))**2.)
        sig_31 =  (0.65*ppicm+0.35*ppi2cm)/2./kcm/
     &           ((W2 - mass(4)**2.)**2. + (mass(4)*width(4))**2.)
        if(sf.EQ.2) then
          width(5)=intwidth(5)*
     &     (xval(47)*petacm/petacmr(5)+(1.-xval(5))*ppi2cm/ppi2cmr(5))

          sig_32 =  (xval(47)*petacm+(1.-xval(47))*ppi2cm)/2./kcm/
     &           ((W2 - mass(5)**2.)**2. + (mass(5)*width(5))**2.)

          width(6)=intwidth(6)*
     &     (xval(48)*petacm/petacmr(5)+(1.-xval(48))*ppi2cm/ppi2cmr(5))

          sig_4 =  (xval(48)*petacm+(1.-xval(48))*ppi2cm)/2./kcm/
     &           ((W2 - mass(6)**2.)**2. + (mass(6)*width(6))**2.)

        endif
        

        sig_del = height(1)*sig_del
        sig_21 = height(2)*sig_21
        sig_22 = height(3)*sig_22
        sig_31 = height(4)*sig_31
        sig_32 = height(5)*sig_32
        sig_4   = height(6)*sig_4

        sig_nr = 0.

        if(sf.EQ.1) then

          do i=1,2  
            sig_nr = sig_nr +(nr_coef(i,1)*(wdif)**(float(2*i+1)/2.))
     &       /(q2+nr_coef(i,2))**
     &       (nr_coef(i,3)+nr_coef(i,4)*q2+xval(44+i)*q2*q2)
          enddo

          sig_nr = sig_nr*xb


        elseif(sf.EQ.2) then
          do i=1,1 
            sig_nr = sig_nr +nr_coef(i,1)*wdif**(float(i)/2.)*
     &       (1.-xth(1))**2.*q2/(1.+nr_coef(i,3)*q2)
     &       /(1.+ q2/nr_coef(i,2))**nr_coef(i,4)/w2

          enddo                           
        endif

        sig_res = sig_del + sig_21 + sig_22 + sig_31 + sig_32 + sig_4
 
        sig_res = sig_res + sig_4L

c        if(sf.EQ.2) then
c          sig_nr = sig_nr*q2/(1.+xval(49)*q2)
c        endif

        sig = sig_res + sig_nr

c        sig = sig_res  

        if(w2.LE.(mp+mpi)**2.OR.sig.LT.0) sig = 0.d0
          
        if(L.EQ.1) sigtemp = sig  

      enddo
       
      q2 = q2temp

      if(lowq2) then
        if(sf.EQ.1) then
          slope = (sigtemp - sig)/0.1
          sig = sigtemp + slope*(q2low-q2)
        else
          slope = sig/q2low
          sig = sig - slope*(q2low-q2)     
        endif
      endif
c      if(lowq2) write(6,*) q2, sig,sigtemp,slope

c      if(sf.eq.1.AND.q2.GT.5) write(6,1000) sig,sig_res,sig_nr 

 1000  format(9f12.5)

      RETURN 
      END 



