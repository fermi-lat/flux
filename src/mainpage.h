// $Header: /nfs/slac/g/glast/ground/cvs/FluxSvc/src/mainpage.h,v 1.11 2003/05/15 20:34:42 burnett Exp $
// Mainpage for doxygen

/*! \mainpage package flux

   \authors Toby Burnett, Sean Robinson, Theodore Hierath, and others.

 \section intro Introduction
  This package implements source definitions, and a procedure to extend them.
  <br>


 <hr>
  @section Basic_XML_Sources Sources
  This is a limited selection of sources. See the contents of the source library file for a complete list.
  @param default
  0.1 GeV gamma-rays coming from the vertical local direction.  Used for default tests.
  @param albedo_gamma
  Source that represents the Earth horizon albedo 
  @param albedo_electronpositron
  Source that represents the splash and re-entrant albedo electrons and positrons
  @param diffuse
  diffuse extraglactic from 10 MeV: from APJ 494:523
  @param diffuse-100mev
  Diffuse extraglactic from 100 MeV
  @param crab-galactic
  the Crab, pulsed portion, with pointed observation, for photons above 100 MeV,
  based on Nolan APJ 409:697
  @param electron
  galactic electron spectrum
  @param normal_gamma
  E^-1 spectrum from 18 MeV to 18 GeV and normal incidence
  @param cosmic_muons
  special source that mimics non-interacting cosmics

 <hr>
 The current complete list is:
 @verbatim
 albedo_electronpositron
albedo_electronpositron_total
albedo_electronpositronavg_total
albedo_electronpositronavghi
albedo_electronpositronavglow
albedo_electronpositronhi
albedo_electronpositronlow
albedo_gamma
albedo_proton_avg
albedo_proton_max
all_gamma
backgndavgpdr
backgndmaxpdr
backgndmix
chime
chimeavg
chimemax
chimemin
cosmic_muons
crab-galactic
crab-pulsed-pointed
cremeavg
default
diffuse
diffuse-100mev
electron
electronavg
electronmax
galcenter
gamma_100_gev_normal
gamma_100_gev_uniform
gamma_100_mev_uniform
gamma_10_gev_uniform
gamma_1_gev_30deg
gamma_1_gev_60deg
gamma_1_gev_normal
gamma_1_gev_uniform
normal_gamma
surface_muons
timetick
timetick30s
vela
vertical_muons
vertical_surface_muons
@endverbatim

This is extracted from 
@verbinclude "source_library.xml"

  <br>
  <h2> Defining an external source </h2>
    See the interface definition ISpectrumFactory for information on how to link code external to this package.


  <hr>
  \section notes release notes
  release.notes
  \section requirements requirements
  \include  requirements
  <hr>
   <hr>
  \todo Complete and recalibrate the CompositeDiffuse structure

*/

