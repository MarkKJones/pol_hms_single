      subroutine dc_effcorr(xfp,dc_eff)
c      program cer_effcorr

c      real function cer_effcorr(hsdelta)

      implicit none

      integer*4 bin
      real*8 eff(10),dc_eff,xfp

ccc   Calculate the bin # for the event  ccc
     
ccc THis should be right to left of the plot.
c      data eff/0.51,0.70,0.81,0.80,0.94,0.85,0.84,0.84,0.90,0.85/ !from fit on above for P run 72436
c      data eff/0.78,0.86,0.90,0.92,0.92,0.93,0.94,0.94,0.94,0.94/ !from fit on above for P run 72505
c      data eff/0.71,0.79,0.84,0.86,0.88,0.90,0.92,0.92,0.93,0.93/ !from fit on above for P run 72790
      data eff/0.93,0.93,0.92,0.92,0.90,0.88,0.86,0.84,0.79,0.71/ !from fit on above for P run 72790

c      data eff/0.85,0.90,0.92,0.94,0.94,0.94,0.94,0.94,0.94,0.93/ !from fit on above for P run 72835

c      data eff/0.70,0.78,0.84,0.87,0.88,0.90,0.92,0.93,0.93,0.92/ ! for E run 72782
c for the run range 72782-72862, the dc maps from 72835

c      bin = 2*int(hsdelta + 8.)+1

c      if(abs(hsdelta).LE.8.) then
c         cer_eff = eff(bin)
c      else
c         cer_eff = 1.0
c      endif


*******************
      if((xfp.lt.37.5).and.(xfp.ge.30.0)) bin=1   ! 13th scin
      if((xfp.lt.30.0).and.(xfp.ge.22.5)) bin=2    ! 12th scin
      if((xfp.lt.22.5).and.(xfp.ge.15.0)) bin=3
      if((xfp.lt.15.0).and.(xfp.ge.07.5)) bin=4
      if((xfp.lt.07.5).and.(xfp.ge.00.0)) bin=5
      if((xfp.lt.00.0).and.(xfp.ge.-07.5)) bin=6
      if((xfp.lt.-07.5).and.(xfp.ge.-15.0)) bin=7
      if((xfp.lt.-15.0).and.(xfp.ge.-22.5)) bin=8
      if((xfp.lt.-22.5).and.(xfp.ge.-30.0)) bin=9
      if((xfp.lt.-30.0).and.(xfp.gt.-37.5)) bin=10  ! 4th scin

      if(abs(xfp).lt.37.5) then
         dc_eff = eff(bin)
      else
         dc_eff = 1.0
      endif

      return
      end
