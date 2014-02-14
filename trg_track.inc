      integer max_z,max_r       !used for array dimensions
      parameter (max_z=700,max_r=700)

      integer num_z(max_r),num_r !actual number of vals from field map

      real*8 rr(max_r)          !array of radial coordinates
      real*8 zz(max_z,max_r)    !array of axial coords, for each radial value
      REAL*8 B_map_z(max_z,max_r),B_map_r(max_z,max_r)
      REAL*8 B_theta_e,B_stheta_e,B_ctheta_e,B_phi_e,B_sphi_e,B_cphi_e 
      REAL*8 B_theta_p,B_stheta_p,B_ctheta_p,B_phi_p,B_sphi_p,B_cphi_p 

      COMMON  /trgFieldStrength/ B_map_z,B_map_r,zz,rr,num_z,num_r
      COMMON  /trgFieldAngles_e/ B_theta_e,B_stheta_e,B_ctheta_e,
     >                           B_phi_e,  B_sphi_e,  B_cphi_e 
      COMMON  /trgFieldAngles_p/ B_theta_p,B_stheta_p,B_ctheta_p,
     >                           B_phi_p,  B_sphi_p,  B_cphi_p 
