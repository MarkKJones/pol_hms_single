pol_hms_single
==============

RSS Hall C single arm Monte Carlo 

The CERN library definition is hardwired into Makefile, so
need to modify.

To compile type "make"

Main code is mc_hms_single.f

=Sub Directories=
*hms  : HMS subroutines
*infiles : input files
*tarfiles : info that describes the target
*outfiles : the output file 
*outfiles/hbook : hbook file

=Info on infiles=
* Mostly self explanatory
* dp/p down and up should be at least -12.0 and 12.0 
* theta down and up is dy/dz , horizontal angle relative to central ray keep at least -40,40
* phi down and up is the dy/dz, vertical angle relative to central ray , this depends on the bending of the particle
* Remember to keep the Dp/p,theta,phi and ztgt reconstruction cut larger than thrown. 
* Keep QFS internal radiation flag =0 , nut used.
* packing fraction=packfrac . If target is not helium filled NH or ND3 
then keep at 0.5. If target is helium filled NH3 or ND3, then 
it is assumed that N and H are the first two materials and helium is
the third. The density of the first 2 materials in the target file is multipled  by packfrac/0.5 and the helium by (1-packfrac/0.5).  