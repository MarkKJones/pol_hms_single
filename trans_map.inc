      INTEGER NMAPS,NLINES 
      PARAMETER (NMAPS  = 20)
      PARAMETER (NLINES = 10000)
 
      INTEGER first(NMAPS),last(NMAPS),order
      INTEGER e1(NLINES),e2(NLINES),e3(NLINES) 
      INTEGER e4(NLINES),e5(NLINES),e6(NLINES)  
      REAL*8    c1(NLINES),c2(NLINES),c3(NLINES) 
      REAL*8    c4(NLINES),c5(NLINES) 
      
      COMMON  /hmsTrack/first,last,
     >                  e1,e2,e3,e4,e5,e6,c1,c2,c3,c4,c5,order

      REAL*8    delta,mep0
      COMMON  /hmsTrackLabCOSY/delta,mep0
