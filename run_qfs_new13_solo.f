      program plot_qfs
c
      implicit none
      real*8 z,a,e,th,w,delnu,elflag
      logical*4 flag
      real*8 sigh,sign,sighe
      real*8 scalein,q2,enu,mp
      real*8 f,tarlen,caplen,rho_nh3
      real*8 trho_h,trho_n,trho_he
      real*8 df_born,df_born2,scale_in
      integer i
c
      e=4730.0
      delnu=14.0
      flag=.false.
      elflag=0.d0
      scale_in = 1.0d0
      mp=.938*1000.
      q2=3.0*1000.*1000.
*
      f = 0.68
      tarlen=3.
      caplen=1.
      rho_nh3=0.867
      trho_h=tarlen*f*rho_nh3*3./17.
      trho_n=tarlen*f*rho_nh3*14./17./14.
      trho_he=(tarlen*(1-f)+caplen)*0.1450/4.
      sigh=0.0d0
*
      do i=1,50
         w=(1.1+(i-1)*.025)*1000.
         enu=(w*w+q2-mp*mp)/2./mp
         th=2*asin(sqrt(q2/4./e/(e-enu)))
         th=th*180.0/3.14159
         z=1.d0
         a=1.d0
        call qfs(z,a,e,th,enu,delnu,flag,sigh,elflag
     >,scale_in)
         z=2.d0
         a=4.d0
        call qfs(z,a,e,th,enu,delnu,flag,sighe,elflag
     >,scale_in)
         z=7.d0
         a=14.d0
        call qfs(z,a,e,th,enu,delnu,flag,sign,elflag
     >,scale_in)
      rho_nh3=0.867
      trho_h=tarlen*f*rho_nh3*3./17.
      trho_n=tarlen*f*rho_nh3*14./17./14.
      trho_he=(tarlen*(1-f)+caplen)*0.1450/4.
         df_born=(trho_h*sigh)/(trho_h*sigh+trho_n*sign+trho_he*sighe)
      rho_nh3=0.867
      trho_h=tarlen*f*rho_nh3*3
      trho_n=tarlen*f*rho_nh3
      trho_he=(tarlen*(1-f)+caplen)
         df_born2=(trho_h*sigh)/(trho_h*sigh+trho_n*sign+trho_he*sighe)
         write(*,*) q2,w,e,e-enu,th,df_born,df_born2,sigh
         enddo
c
      end
