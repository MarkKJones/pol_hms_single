! STRUCT_HMS.INC

!make fortran 77 compatible by removing structure declaration
!use equivalence statement to link hstop array to indiv hstop_xxx vars

      integer hSTOP(23)
      integer hSTOP_trials,hSTOP_slit_hor,hSTOP_slit_vert
     +  ,hSTOP_slit_oct,hSTOP_Q1_in,hSTOP_Q1_mid,hSTOP_Q1_out
     +  ,hSTOP_Q2_in,hSTOP_Q2_mid,hSTOP_Q2_out,hSTOP_Q3_in
     +  ,hSTOP_Q3_mid,hSTOP_Q3_out,hSTOP_D1_in,hSTOP_D1_out,hSTOP_hut
     +  ,hSTOP_dc1,hSTOP_dc2,hSTOP_scin,hSTOP_cal,hSTOP_successes
     +  ,nfail_track_to_tgt,nfail_track_from_tgt

      equivalence
     +  (hSTOP( 1),hSTOP_trials),
     +  (hSTOP( 2),hSTOP_slit_hor),
     +  (hSTOP( 3),hSTOP_slit_vert),
     +  (hSTOP( 4),hSTOP_slit_oct),
     +  (hSTOP( 5),hSTOP_Q1_in),
     +  (hSTOP( 6),hSTOP_Q1_mid),
     +  (hSTOP( 7),hSTOP_Q1_out),
     +  (hSTOP( 8),hSTOP_Q2_in),
     +  (hSTOP( 9),hSTOP_Q2_mid),
     +  (hSTOP(10),hSTOP_Q2_out),
     +  (hSTOP(11),hSTOP_Q3_in),
     +  (hSTOP(12),hSTOP_Q3_mid),
     +  (hSTOP(13),hSTOP_Q3_out),
     +  (hSTOP(14),hSTOP_D1_in),
     +  (hSTOP(15),hSTOP_D1_out),
     +  (hSTOP(16),hSTOP_hut),
     +  (hSTOP(17),hSTOP_dc1),
     +  (hSTOP(18),hSTOP_dc2),
     +  (hSTOP(19),hSTOP_scin),
     +  (hSTOP(20),hSTOP_cal),
     +  (hSTOP(21),hSTOP_successes),
     +  (hSTOP(22),nfail_track_to_tgt),
     +  (hSTOP(23),nfail_track_from_tgt)

      common /spectrometer1/ hSTOP
