*------------------------------------------------------------------------
*
*       GEN_TRACK  Gen Tracking Routines 
*      -=========-
* 
*	Backward Tracking of electrons in the Jlab HMS hall 
*       C spectrometer including the Gen's target magnetic field
*
*       Note: - the GEN routines use a lab (HMS) coord. system
*               and a BEAM coord. system, both right handed with
*                 x : pointing downwards
*                 y : perpendicular to x,z, 
*                     pointing to the left (if seen in z-direction)
*                 z : HMS or BEAM axis, pointing from the target to the 
*                     focal plane or downstream, respectively
*
*             - all lengths (x,y,z,l,...) are measured in [m]
*             - angles are measured either as dx/dz,dy/dz (HMS coords.)
*               or counter clock wise in [deg]     
*             - the momentum is measured in delta (relative momentum
*               deviation  = 1-p0/pHMS)
*
*       Supplies:
*         genInitRecon (rmap,fmap,thetaB,cmap,order,thetaE,p0)
*           load the target field map and the HMS forward (cosy) and 
*           reconstruction maps 
*         genRecon (u,x,y,uT)
*           reconstruction of the target coordinates 
*           (delta, dx/dz, y, dy/dz) including the effects of the
*           target magnetic field and the vertical beam offset
*       
*       Note: - Before calling genRecon the target field map, the forward 
*               and backward maps have to be loaded by a call to 
*               genInitRecon
*
*       Requires: 
*         the hms_track and trg_track modules
*
*       written by Markus Muehlbauer for the GEN Experiment
*------------------------------------------------------------------------

      SUBROUTINE genInitRecon (recon,field,thetaB,cosy,order,thetaE,p0) 
      IMPLICIT  NONE
      CHARACTER recon*(*),field*(*),cosy*(*)
      REAL      thetaE,thetaB,p0
      INTEGER   order
      
      include 'gen_constants.par'
      
* --  load the target field map and the HMS forward (cosy) map 
*
*     Parameter:
*       recon      I : name of the reconstruction map
*       field      I : name of the target field map
*       thetaB     I : angle of the magnetic field [deg]
*       cosy       I : name of the forward map
*       order      I : maximal order to take into account
*       thetaE     I : spectrometer angle [deg]
*       p0         I : momentum the spectrometer is set to

      REAL    theta,ctheta,stheta,p 
      COMMON /genParameter/theta,ctheta,stheta,p 

      theta  = thetaE
      ctheta = COS(thetaE*degree)
      stheta = SIN(thetaE*degree)
      p      = p0
      
      CALL hmsInitRecon   (recon)
      CALL trgInit        (field,thetaB-thetaE,0.)
      CALL hmsInitForward (cosy,order,.FALSE.,p0) 
       
      RETURN
      END
      
*------------------------------------------------------------------------
      
      SUBROUTINE genRecon (u,x,y,uT,ok)
      IMPLICIT NONE
      REAL     u(4),x,y,uT(6)
      LOGICAL  ok
            
      include 'gen_constants.par'

* --  performs the reconstruction of the target coordinates 
*     (delta, dx/dz, y, dy/dz) including the effects of the
*     target magnetic field and the vertical beam offset
*     
*     Parameter:
*       u      I : focal plane coordinates  
*                    u(1,2)  : x [m], dx/dz = out of plane coords. (downwards) 
*                    u(3,4)  : y [m], dy/dz = inplane coords. (perp. on x,z)
*       x      I : vert. beam offset [m] (out of plane coord.; downwards)
*       y      I : hori. beam offsey [m] (inplane coord.; perp on x-beam, z-beam)
*       uT     O : target coordinates
*                    uT(1,2) : x [m], dx/dz = out of plane coord. (downwards) 
*                    uT(3,4) : y [m], dy/dz = inplane coord. (perp. on x,z)
*                    uT(5)   : z [m] = in axis coordinate (towards HMS)  
*                    uT(6)   : delta = relative deviation of the particle 
*                                      momentum from p0
*       ok   IO  : status variable 
*                   - if false no action is taken 
*                   - set to false when no reconstruction is found 

      REAL theta,ctheta,stheta,p
      COMMON /genParameter/theta,ctheta,stheta,p 
     
      REAL    xx,dx,vT(6),vTx(6)
      INTEGER i,n
      
      real REF_VAL
      parameter (REF_VAL=100.)         ! what does this acyually correspond to?
      real OTHER_REF
      parameter (OTHER_REF=30.)         ! what does this acyually correspond to?

      REAL       eps        ! accurracy for x in mm
      PARAMETER (eps = 0.2) ! (one more iteration is performed 
                            !  after the given accuraccy is reached) 
      xx = x
      
      ! find a first approximation for uT
      CALL hmsReconOutOfPlane (u,xx,uT,ok)
      IF (.NOT. ok) RETURN
   	 
      ! drift to a field free region and calculate the velocities
      vT(1) = REF_VAL*(uT(1)+1.*uT(2))
      vT(2) = REF_VAL*(uT(3)+1.*uT(4))
      vT(3) = REF_VAL*1.
      vT(6) = OTHER_REF/SQRT(1+uT(2)**2+uT(4)**2)
      vT(4) = uT(2)*vT(6)
      vT(5) = uT(4)*vT(6) 
      
      ! and track into the magnetic field to the beam plane (perp. to y)
      CALL trgTrackToPlane (vT,-p*(1+uT(6))/MeV,1.,
     >                      0.,-ctheta,stheta,y*REF_VAL,ok)
    
      n  = 0
      dx = 1.
      DO WHILE ((dx .GT. .2) .AND. (n .LT. 10) .AND. ok)
        dx = abs(x*REF_VAL-vT(1))  
        ! track to the z=0 plane to find a correction for the x-offset   
        vTx(1) = REF_VAL*x 
        DO i=2,6
          vTx(i) = vT(i)
        ENDDO
    
        CALL trgTrackToPlane (vT, -p*(1+uT(6))/MeV,1.,0., 0.,1.,0., ok)
        CALL trgTrackToPlane (vTx,-p*(1+uT(6))/MeV,1.,0., 0.,1.,0., ok) 
               
        xx = xx+(vTx(1)-vT(1))*0.01   ! what unit conversion is this???
      
        ! now find a better approximation for uT
        CALL hmsReconOutOfPlane (u,xx,uT,ok)
         
        ! drift to a field free region and calculate the velocities
        vT(1) = REF_VAL*(uT(1)+1.*uT(2))
        vT(2) = REF_VAL*(uT(3)+1.*uT(4))
        vT(3) = REF_VAL*1.
        vT(6) = OTHER_REF/SQRT(1+uT(2)**2+uT(4)**2)
        vT(4) = uT(2)*vT(6)
        vT(5) = uT(4)*vT(6) 
     
        ! and track into the magnetic field to the beam plane (perp. to y)
        CALL trgTrackToPlane (vT,-p*(1+uT(6))/MeV,1.,
     >                        0.,-ctheta,stheta,y*REF_VAL,ok) 
        n = n+1
      ENDDO
  
      IF (dx .GT. .2) ok = .FALSE.
            
      ! calculate the result in HMS coordinates
      uT(1) =  0.01*vT(1)  
      uT(2) = vT(4)/vT(6)
      uT(3) =  0.01*vT(2)
      uT(4) = vT(5)/vT(6)
      uT(5) =  0.01*vT(3)
   
      RETURN
      END
 
      
