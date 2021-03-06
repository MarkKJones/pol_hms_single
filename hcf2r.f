!     subroutine h2model_liang(qq,ww,w1,w2)
      subroutine hcf2r(qq,ww,f2,r)
                       
c     include 'model_liang.cmn'
      implicit none  
      integer i,j,j1,i1,tmp
      real*8 qq,ww,w1,w2,ww1,ww2  !,w11,w22
      real*8 r,f2,nu,mp,mp2,r1,f21,r2,f22
      real*8 q4,q6,xx,fac,rlog,q8 !,frac
      real*8 rp(276), f2p(552),tmp1,tmp2
      real*8 q,qth,dq,y,f2th,lin,binding
      real*8 scb,quad
!	real*8 mp, mp2	
!       real*8 x43(100),q443(100),q743(100)
!	real*8 w
      logical test
      common/model_liang/rp,f2p
      common/model_liang1/test
 
      mp =0.938272d0
      mp2 = mp*mp

  
c      write(*,*) qq,ww
      if (ww.lt.1.15d0) then
         w1 = 0.d0
         w2 = 0.d0        
      else if (ww.gt.3.9d0) then
	return
!         call h2model(qq,ww,w1,w2)
      else
c       write(*,*) test
        if (.not.test) then 
c        write(*,*) 'reading parameters for liang model'
        test = .true.


      data rp/
     > 0.33129d0, 1.32645d0,-0.79900d0,  -3.11559d0,    2.5681d0,  -1.82534d0,
     > 0.28813d0, 0.90973d0,-0.54869d0,  -3.76651d0,    3.4846d0,  -2.02225d0,
     > 0.56382d0, 0.43268d0, 0.08276d0,  23.37011d0,  -25.4063d0,  13.71927d0,
     > 0.67816d0,-0.09742d0, 1.91168d0,  51.09415d0,  -61.0272d0,  17.29566d0,
     > 0.64090d0, 0.51581d0, 0.28952d0, -69.99436d0,  -98.0096d0, 133.50420d0,
     > 0.46245d0,-0.02240d0, 0.31779d0,  -3.26504d0,   77.9481d0,  14.49329d0,
     > 0.39436d0,-0.02232d0, 0.34328d0,  24.69035d0,   50.0479d0,  13.91268d0,
     > 0.39805d0, 0.03919d0, 0.07716d0,   3.75484d0,   10.3324d0,   8.08673d0,
     > 0.38997d0, 0.00900d0, 0.15343d0,  16.36940d0,    4.9903d0,   7.44212d0,
     > 0.29674d0, 0.02783d0, 0.06118d0,  27.36369d0,  -17.4484d0,   4.32732d0,
     > 0.35439d0,-0.01126d0, 0.23150d0,  53.46022d0,  -18.2759d0,   7.18313d0,
     > 0.31987d0, 0.01455d0, 0.23862d0,  81.55396d0,  -74.7982d0,   5.96615d0,
     > 0.31163d0,-0.00699d0, 0.12995d0,  51.23706d0,   18.1217d0,   4.82674d0,
     > 0.32176d0, 0.12831d0, 0.00514d0,  68.39391d0,  -73.3401d0,   3.86676d0,
     > 0.32514d0, 0.20568d0, 0.02398d0,  31.46726d0,  -34.1291d0,   3.42527d0,
     > 0.48734d0, 0.32263d0, 0.07471d0,  42.71339d0,  -49.2712d0,   7.19561d0,
     > 0.38366d0, 0.16378d0, 0.00110d0,  26.58418d0,  -29.0223d0,   3.80760d0,
     > 0.52817d0,-0.37320d0,10.22177d0,  44.87227d0,  -53.3864d0,  15.03388d0,
     > 0.34971d0, 0.24432d0, 0.03942d0,  15.46841d0,  -17.0274d0,   2.87509d0,
     > 0.32309d0, 0.05971d0, 0.23399d0,  81.18824d0,  -85.6957d0,   4.59427d0,
     > 0.48408d0, 0.14606d0,-0.12127d0,  78.59016d0,  -92.2077d0,   2.15488d0,
     > 0.44745d0, 0.08725d0,-0.07341d0,  96.91505d0, -113.3328d0,   1.61255d0,
     > 0.17621d0, 0.00267d0, 0.16446d0, 290.88310d0, -284.5832d0,   4.78794d0,
     > 0.19235d0, 0.03305d0, 0.26250d0,  86.32166d0,  -86.1215d0,   4.45001d0,
     > 0.25966d0, 0.02309d0, 0.06332d0,  54.29687d0,  -49.3171d0,   3.33319d0,
     > 0.35243d0, 0.02105d0,-0.18868d0, 110.69520d0, -226.4427d0,   4.95000d0,
     > 0.36534d0,-0.21314d0, 2.81021d0, -59.05300d0,  -37.7424d0,  31.42299d0,
     > 0.28897d0,-5.27934d0,62.75448d0,  15.15152d0,    5.7736d0,  29.82406d0,
     > 0.31192d0,-0.03546d0, 0.09286d0,  19.06862d0,  -33.4685d0,   1.56198d0,
     > 0.35713d0,-0.20097d0, 2.05940d0, 436.08360d0, -517.4285d0,  14.75464d0,
     > 0.32423d0, 0.05549d0, 0.04984d0,  74.70284d0,  -85.5771d0,   3.61358d0,
     > 0.26321d0, 0.02532d0, 0.04165d0, 291.81850d0, -330.3480d0,   3.42817d0,
     > 0.19235d0,-0.00521d0, 0.30864d0, 262.23130d0, -271.3294d0,   4.52254d0,
     > 0.17858d0,-0.00018d0, 0.28614d0, 195.28650d0, -195.6441d0,   3.92928d0,
     > 0.14180d0, 0.00801d0, 0.21545d0, 131.05210d0, -121.0251d0,   3.20372d0,
     > 0.19810d0, 0.00240d0, 0.29005d0, 122.32390d0, -121.0733d0,   3.60232d0,
     > 0.24832d0,-0.00999d0, 0.46397d0,  81.22332d0,  -85.5181d0,   3.40233d0,
     > 0.31740d0, 0.12344d0, 0.07026d0,  53.01398d0,  -64.9467d0,   3.08821d0,
     > 0.32739d0,-0.06064d0, 0.79821d0, 104.36790d0, -124.3242d0,   4.74774d0,
     > 0.34597d0, 0.11624d0,-0.04088d0,  50.29695d0,  -64.8436d0,   2.35262d0,
     > 0.33840d0, 0.11035d0, 0.03128d0,  39.27810d0,  -51.3317d0,   2.38903d0,
     > 0.38341d0, 0.12239d0,-0.12886d0,  59.12659d0,  -80.0135d0,   1.94849d0,
     > 0.34759d0,-0.06730d0, 0.99427d0,  87.77917d0, -116.0158d0,   4.31371d0,
     > 0.35035d0, 0.13825d0, 0.05372d0,  37.00152d0,  -51.0832d0,   2.35058d0,
     > 0.31877d0, 0.09654d0, 0.12551d0,  70.29255d0,  -93.9217d0,   2.72640d0,
     > 0.32985d0, 0.25111d0,-0.24644d0,  11.54706d0,  -17.2718d0,   1.00506d0/

      data f2p/
     >   0.00270d0,  -0.02563d0,  0.06438d0,   0.13288d0,  -0.18442d0,  59.57852d0,
     > 126.68970d0,-156.52300d0,259.6862d0,   45.71608d0,  51.94866d0,-337.30260d0,
     >  -0.01266d0,   0.53825d0, -1.14310d0,  -1.10755d0,   2.67523d0, -96.68052d0,
     > -65.71029d0, 144.62170d0,-350.7377d0,  61.76398d0, 116.55460d0, 172.52990d0,
     >  -0.00820d0,   0.30865d0, -0.10750d0,  -0.79825d0,   1.01081d0, -37.89112d0,
     > -99.42540d0, 120.79090d0,-161.3490d0, -27.73087d0, -37.25240d0, 223.31690d0,
     >  -0.01626d0,   0.68263d0, -2.08330d0,   0.93416d0,   7.66115d0, -20.34949d0,
     > -89.79055d0, 100.32260d0, -89.9365d0, -81.08472d0,  33.33122d0, 139.55320d0,
     >   0.10085d0,  -0.39180d0,  4.39375d0,   2.96299d0,  -9.43805d0,   5.53277d0,
     > -49.30489d0,  43.73740d0,  12.3122d0,-123.36040d0, 117.44570d0,  -2.55460d0,
     >   0.16801d0,   1.66346d0, -3.58292d0,  -1.43295d0,   6.54949d0,   9.59234d0,
     > -15.43190d0,   7.56539d0,  32.3811d0, -79.52365d0,  79.85619d0, -30.67446d0,
     >   0.24036d0,   0.66839d0, -1.02483d0,   2.73454d0,  -3.12791d0,   4.67220d0,
     >  -1.94161d0,  -1.02171d0,  14.8272d0, -21.21147d0,  20.23464d0, -11.97409d0,
     >   0.14395d0,  -0.24303d0,  0.47600d0,   5.20955d0,   4.07727d0,  -6.50083d0,
     > -59.81225d0,  63.16642d0, -30.5782d0, -70.20249d0,  50.86131d0,  55.80522d0,
     >   0.19595d0,   0.49078d0,  0.09058d0,   2.08639d0,   0.61890d0,   3.84467d0,
     >  30.61177d0, -30.33699d0,  17.4414d0,  32.96692d0, -20.06578d0, -29.88399d0,
     >   0.16363d0,   0.21846d0,  0.20264d0,   2.53441d0,  -0.76028d0,   3.90782d0,
     >  36.58707d0, -35.68809d0,  18.2261d0,  40.25367d0, -24.84555d0, -33.50865d0,
     >   0.11976d0,  -0.45730d0,  1.08627d0,   2.63121d0,  -0.58440d0,  -5.72117d0,
     > -53.80934d0,  57.73899d0, -25.7978d0, -52.77444d0,  37.38223d0,  48.12931d0,
     >   0.03129d0,   0.08877d0,  0.31201d0,   0.52214d0,   0.56099d0,   9.45694d0,
     >  80.42830d0, -77.88046d0,  41.7567d0,  69.77966d0, -24.70381d0, -89.21936d0,
     >   0.04149d0,   0.05422d0,  0.14645d0,   0.46049d0,   0.70058d0,   7.52857d0,
     >  73.36507d0, -69.99093d0,  33.8845d0,  68.44820d0, -34.26609d0, -69.93432d0,
     >   0.03743d0,   0.05538d0,  0.06989d0,   0.42875d0,   0.74432d0,   9.91687d0,
     >  99.10915d0, -95.21824d0,  44.3231d0,  88.44396d0, -44.56128d0, -92.21440d0,
     >   0.04214d0,   0.07361d0, -0.05129d0,   0.40464d0,   1.03509d0,   8.78332d0,
     >  90.67270d0, -86.53600d0,  39.0777d0,  79.55918d0, -42.40906d0, -79.76051d0,
     >   0.02977d0,   0.24989d0, -0.36480d0,   0.43282d0,   1.66379d0,  10.58702d0,
     >  99.35179d0, -95.66632d0,  45.5028d0,  74.12631d0, -26.76438d0, -97.86623d0,
     >   0.04183d0,   0.17483d0, -0.23191d0,   0.45452d0,   1.41034d0,   8.17936d0,
     >  74.28743d0, -69.96535d0,  34.5675d0,  50.35936d0, -11.35428d0, -76.36961d0,
     >   0.05173d0,   0.04448d0,  0.13424d0,   0.38647d0,   0.54959d0,   6.16662d0,
     >  51.08549d0, -46.36961d0,  25.3451d0,  28.23540d0,   3.70411d0, -58.01252d0,
     >   0.08047d0,   0.07636d0, -0.05934d0,   0.34274d0,   0.77572d0,   4.27257d0,
     >  38.90837d0, -35.54622d0,  17.7386d0,  24.76801d0,  -9.17321d0, -33.32258d0,
     >   0.08916d0,   0.08720d0, -0.04772d0,   0.36316d0,   0.77081d0,   4.10228d0,
     >  40.26231d0, -37.09061d0,  17.1693d0,  26.59951d0, -14.07698d0, -29.96852d0,
     >   0.02688d0,   0.20865d0,  0.24902d0,   0.19521d0,   0.37927d0,  -0.57453d0,
     >  -9.00651d0,  14.38374d0,  -2.4497d0,  -6.51968d0,   9.73428d0,   3.89960d0,
     >   0.04139d0,   0.15301d0,  0.32248d0,   0.36027d0,   0.17856d0,   3.32167d0,
     >  28.88752d0, -23.90850d0,  13.4160d0,  12.82948d0,   3.90580d0, -29.60664d0,
     >   0.05169d0,  -0.19792d0,  0.70477d0,   0.98245d0,   0.19970d0,   1.07366d0,
     >  -0.96632d0,  13.46862d0,   3.5916d0, -10.17772d0,  41.48280d0, -29.53078d0,
     >   0.04810d0,  -0.14599d0,  0.44319d0,   0.86999d0,   0.59196d0,   2.37207d0,
     >  16.49698d0,  -3.58537d0,   9.1814d0,   2.29473d0,  31.25257d0, -39.15542d0,
     >   0.04731d0,  -0.14069d0,  0.33899d0,   0.74002d0,   0.62827d0,  -0.83045d0,
     > -17.17104d0,  32.13834d0,  -3.8414d0, -13.50958d0,  37.48258d0, -11.74448d0,
     >   0.05221d0,  -0.11177d0,  0.27687d0,   0.68277d0,   0.62925d0,   2.10930d0,
     >  13.86561d0,  -1.66573d0,   7.9977d0,   0.35896d0,  27.99104d0, -32.29613d0,
     >   0.06442d0,  -0.07257d0,  0.11952d0,   0.50396d0,   0.53317d0,   4.29581d0,
     >  44.07500d0, -35.61928d0,  17.2312d0,  19.12265d0,   0.86142d0, -37.08731d0,
     >   0.07351d0,  -0.10911d0,  0.12156d0,   0.46602d0,   0.48331d0,   5.17804d0,
     >  58.31429d0, -50.40657d0,  21.0203d0,  27.73335d0, -10.71966d0, -39.46624d0,
     >   0.06943d0,  -0.07336d0,  0.15002d0,   0.46675d0,   0.51645d0,   3.31844d0,
     >  37.30048d0, -28.68561d0,  13.4842d0,  17.16954d0,  -2.77701d0, -26.91567d0,
     >   0.05361d0,  -0.02633d0,  0.19376d0,   0.44691d0,   0.44814d0,   2.88220d0,
     >  29.95785d0, -19.93872d0,  11.4738d0,  10.83239d0,   7.47675d0, -28.12408d0,
     >   0.03795d0,  -0.04353d0,  0.30496d0,   0.53973d0,   0.44511d0,   1.82131d0,
     >  14.62905d0,  -0.35599d0,   6.9752d0,   1.48548d0,  24.99437d0, -29.88984d0,
     >   0.04790d0,   0.05632d0,  0.20660d0,   0.39304d0,   0.45673d0,   4.19462d0,
     >  47.41709d0, -39.83233d0,  16.6487d0,  17.18184d0,  -2.52659d0, -32.65720d0,
     >   0.05140d0,   0.11201d0,  0.19756d0,   0.26992d0,   0.30311d0,   3.20459d0,
     >  37.47848d0, -31.60155d0,  12.7851d0,  13.97016d0,  -8.18266d0, -18.63125d0,
     >   0.04703d0,   0.17203d0,  0.21807d0,   0.23841d0,   0.26405d0,   1.17976d0,
     >  10.84169d0,  -4.87255d0,   4.5850d0,   2.60271d0,  -1.08836d0,  -2.73421d0,
     >   0.09651d0,  -0.07097d0, -0.01140d0,   0.41291d0,   0.67944d0,   8.63177d0,
     > 116.48230d0,-115.22170d0,  34.6775d0,  50.71781d0, -62.71606d0, -31.72991d0,
     >   0.05407d0,  -0.01853d0,  0.18484d0,   0.39172d0,   0.39883d0,   4.54709d0,
     >  57.66357d0, -47.37784d0,  18.0009d0,  22.01861d0, -16.19442d0, -25.28908d0,
     >   0.05646d0,  -0.02226d0,  0.19051d0,   0.38210d0,   0.35937d0,   2.35732d0,
     >  31.71869d0, -21.30330d0,   9.5296d0,  14.48209d0, -14.30969d0,  -7.90519d0,
     >   0.04477d0,   0.00480d0,  0.15076d0,   0.29189d0,   0.31529d0,   3.12081d0,
     >  41.81401d0, -31.97399d0,  12.4581d0,  15.71344d0, -15.68024d0, -13.09181d0,
     >   0.03143d0,   0.06733d0,  0.18032d0,   0.27882d0,   0.27721d0,   1.09240d0,
     >  15.14939d0,  -7.18577d0,   4.4895d0,   7.48708d0, -11.64642d0,   1.23596d0,
     >   0.09641d0,  -0.07037d0, -0.04597d0,   0.26390d0,   0.48362d0,   3.68185d0,
     >  61.90438d0, -66.11752d0,  15.0943d0,  29.88299d0, -76.76054d0,  30.73642d0,
     >   0.05348d0,  -0.06567d0,  0.26431d0,   0.51672d0,  -0.39508d0,  -0.17103d0,
     >   0.38781d0,   9.77783d0,  -0.2076d0,   3.74785d0, -13.89592d0,  15.93324d0,
     >   0.01978d0,   0.02381d0,  0.16282d0,   0.27720d0,   0.28227d0,   1.10671d0,
     >  19.57988d0,  -8.77431d0,   4.6923d0,  10.90089d0, -21.92500d0,   4.82789d0,
     >   0.02816d0,   0.08035d0,  0.15177d0,   0.19596d0,   0.20729d0,   0.16046d0,
     >   4.85393d0,   1.60025d0,   0.9924d0,   6.34038d0, -23.91168d0,  19.81668d0,
     >   0.00810d0,   0.06904d0,  0.19159d0,   0.22518d0,   0.16841d0,   0.22909d0,
     >   4.46808d0,  -1.67027d0,   1.1692d0,   4.17188d0, -22.18001d0,  13.56817d0,
     >  -0.00459d0,   0.08266d0,  0.20067d0,   0.22163d0,   0.15112d0,   0.67785d0,
     >  10.03290d0, -23.97605d0,   2.8133d0,   1.49803d0, -32.05834d0,  14.49983d0,
     >  -0.00869d0,   0.07988d0,  0.20642d0,   0.23439d0,   0.16011d0,   1.58356d0,
     >  23.98258d0, -56.34460d0,   6.2639d0,  -0.04471d0, -43.27986d0,  13.87953d0/
        endif

c      if (ww .ge. 3.9) then
c           ww = 3.875
c      endif

       if (ww.lt.1.2d0) then
        ww1 = 1.175d0
c       ww1 = 1.183d0
        tmp = 0
       else if (ww.lt.3.0d0) then  
        tmp =int((ww-1.15d0)/0.05d0) 
        ww1 = 1.175d0+tmp*0.05d0
       else
        tmp = int((ww-3.0d0)/0.1d0)
        ww1 = 3.05d0+tmp*0.1d0 
        tmp = tmp+37
       endif
        j = (tmp+1)*6-5
        j1 = (tmp+1)*12-11
        q4 = qq*qq
        q6 = qq*qq*qq
        q8 = q4*q4

      tmp1 = 12.d0*(qq/(1.d0+qq))
      tmp2 = log(qq/.04d0)

      xx = qq/(qq+ww1-mp2)
      fac = 1+tmp1*(.125d0**2/(xx**2+.125d0**2))
      rlog = fac/tmp2

      r1 = rp(j)*rlog
     &+(rp(j+1)/qq+rp(j+2)/(q4+.3d0**2))
     &*(1.d0+rp(j+3)*xx+rp(j+4)*(xx**2))
     &*(xx**rp(j+5)) 

      binding = 1.d0
      q = log(qq)                                                       
      qth = .2d0+3.2d0*xx                                                      
      dq = q-qth
      y = 1.d0-xx                                                         
      f2th = f2p(j1)*y**1+f2p(j1+1)*y**2+f2p(j1+2)*y**3+          
     >       f2p(j1+3)*y**4+f2p(j1+4)*y**5                           
      quad = (f2p(j1+5)+f2p(j1+6)*xx+f2p(j1+7)*xx**2)*dq**2      
      lin  = (f2p(j1+8)+f2p(j1+9)*xx+f2p(j1+10)*xx**2
     >+f2p(j1+11)*xx**3)*dq               
      if (q.gt.qth) quad = 0.d0                                              
      scb  = (1.d0+lin+quad)                                                 
      f21  = f2th*scb*binding  
c      write(*,*) f2th,scb,binding

c      if (ww.le.1.2 .or. ww.ge.3.85
c     >.or. abs(ww-ww1).le.0.001) then
      if (ww.le.1.2d0 .or. abs(ww-ww1).le.0.001d0) then
          r = r1
          f2 = f21
      else
          if (ww.gt.ww1) then
             i = j+6
             i1 = j1+12
             if (ww1 .lt. 1.2d0) then
               ww2 = 1.225d0
             else if (ww1 .lt. 2.975d0) then            
               ww2 = ww1+0.05d0
             else if (abs(ww1- 2.975d0) .lt. 0.01d0) then
               ww2 = 3.05d0
             else if (abs(ww1-3.85d0) .lt. 0.01d0) then
               ww2 = 3.75d0
               i = j - 6
               i1 = j1 - 12
             else 
               ww2 = ww1+0.1d0
             endif  
          else 
             i = j-6
             i1 = j1-12
             if (ww1 .lt. 1.25d0) then
               ww2 = 1.175d0
c               ww2 = 1.183d0
             else if (ww1 .lt. 3.05d0) then            
               ww2 = ww1-0.05d0
             else if (abs(ww1-3.05d0) .lt. 0.01d0) then
               ww2 = 2.975d0
             else 
               ww2 = ww1-0.1d0
             endif  
          endif    
      xx = qq/(qq+ww2-mp2)
      fac = 1.d0+tmp1*(.125d0**2/(xx**2+.125d0**2))
      rlog = fac/tmp2

      r2 = rp(i)*rlog
     &+(rp(i+1)/qq+rp(i+2)/(q4+.3d0**2))
     &*(1.d0+rp(i+3)*xx+rp(i+4)*(xx**2))
     &*(xx**rp(i+5))

      q = log(qq)                                                       
      qth = .2d0+3.2d0*xx                                               
      dq = q-qth
      y = 1.d0-xx                                               
      f2th = f2p(i1)*y**1+f2p(i1+1)*y**2+f2p(i1+2)*y**3+          
     >       f2p(i1+3)*y**4+f2p(i1+4)*y**5                                 
      quad = (f2p(i1+5)+f2p(i1+6)*xx+f2p(i1+7)*xx**2)*dq**2          
      lin  = (f2p(i1+8)+f2p(i1+9)*xx+f2p(i1+10)*xx**2
     >+f2p(i1+11)*xx**3)*dq               
      if (q.gt.qth) quad = 0.d0                                              
      scb  = (1.d0+lin+quad)                                                 
      f22  = f2th*scb*binding
c      write(*,*) f2th,scb,binding

       r = r1+(ww-ww1)*(r2-r1)/(ww2-ww1)
       f2 = f21+(ww-ww1)*(f22-f21)/(ww2-ww1)
      endif   

        if (f2.lt.0.d0) then
!           f2 = 0.001d0
           f2 = 0.d0
        elseif (r.lt.0.d0) then
!           r = 0.001d0
           r = 0.d0
        else
        endif  
 
        nu = (qq+ww-mp2)/2.d0/mp
        w2 = f2/nu 
        w1 = w2*(1.0d0 + nu*nu/qq)/(1.0d0 + r)

!        write(*,*) qq,w,ww1,ww2,r,f2,w1,w2        
!	write(logun1,697) w,qq,r,f2,w1,w2
!697	format(6g13.4)
c        if (qq.ge.4.0) then
c           call h2model(qq,ww,w11,w22)
c           frac = (qq - 4.0)/(5.0 - 4.0)
c           w1 = w11*frac + w1*(1.0 - frac)
c           w2 = w22*frac + w2*(1.0 - frac)
c         endif

       endif     

      return
      end                                                               

