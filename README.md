pol_hms_single
==============

RSS Hall C single arm Monte Carlo 

The CERN library definition is hardwired into Makefile, so
need to modify.

To compile type "make"

Main code is mc_hms_single.f

= Running code =

mc_hms_single   infile_name   tarfile_name label

tarfile_name is the target material file in subdirectory : tarfiles
Output file is at outfiles/infile_name_tarfile_name_label.out 
The hbook file is at outfiles/hbook/infile_name_tarfile_name_label.1.hbook 

= Code flow =
* Read in the input file
* Read in target material file
* For each material calculates luminosity = lumin  
*  For each material calculates prob=density*thick
* Starts event loop
* Randomly selects target material weighted by the "prob" factor
* Randomly picks beam raster positin and target position ( depends on selected target material)
* Randomly picks dy/dz and dx/dz
* Calculates the beam energy loss in the target.
* Randoms selects whether the event is elastic scattering or inelastic scattering.
* If elastic calculates scattered energy. If ineasltic randoms picks scattered energy.
* Calculates the scattered particle energy loss and  radiation lengths for radiative code.  
* If radiate flag is true then randomly determines the radiation loss.
* Determines vertex kinematics for cross section calculations. If inelastic 
  and W_vertex < 1.073 then skip event
* If magnet field on then tracks through field
* Tracks through the HMS to see if it makes it to focal plane.
* Call subroutine qfs to get cross section. Code qfs_new13_sub.f is compiled in the Makefile.
**For inelastic uses F1F2IN09 and for quasi-elastic F1F2QE09 and sums the cross sections.
**For hydrogen elastic uses the Bosted fit for the form factors.
* Calculates norm_rate which is put into ntuple variable "weight"


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
* phi down and up is the dx/dz, vertical angle relative to central ray , this depends on the bending of the particle
* Remember to keep the Dp/p,theta,phi and ztgt reconstruction cut larger than thrown. 
* Keep QFS internal radiation flag =0 , not used.
* packing fraction=packfrac . If target is not helium filled NH or ND3 
then keep at 0.5. If target is helium filled NH3 or ND3, then 
it is assumed that N and H are the first two materials and helium is
the third. The density of the first 2 materials in the target file is 
multipled  by packfrac/0.5 and the helium by (1-packfrac/0.5).  

= Target material file =
* List of materials in the target.