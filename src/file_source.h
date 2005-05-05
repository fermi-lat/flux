/**
   @page file_source FileSource 

   This source reads in incident particle properties from a text file
   on an event-by-event basis.  Here is an example input file:

   @verbatim
# test file for FileSource
# columns are
# particle KE(GeV) x(mm) y  z  cos_x cos_y  cos_z
mu-        0.100    0    0    0    1     0      0
mu-        0.300    0    0    0    0     1      0
mu-        0.200    0    0    0    0     0      1
mu- 4.041 -48.48 -38.23 800.00 -0.161 0.227 -0.960
mu- 2.485 280.91 -166.79 800.00 -0.507 0.158 -0.847
mu- 2.292 186.10 105.79 800.00 -0.102 -0.009 -0.995
mu- 9.357 74.38 128.95 800.00 0.086 -0.197 -0.977
mu- 1.777 -45.42 -163.01 800.00 -0.123 0.276 -0.953
mu- 0.642 -71.06 -179.82 800.00 0.229 0.092 -0.969
mu- 2.362 190.36 -100.15 800.00 -0.253 0.297 -0.921
mu- 4.19 -52.47 193.12 800.00 0.113 -0.304 -0.946
mu- 0.607 4.84 -199.59 800.00 0.021 0.160 -0.987
mu- 1.204 -182.71 167.73 800.00 0.091 -0.291 -0.953
mu- 1.783 -120.65 175.54 800.00 0.284 -0.144 -0.948
mu- 2.201 -167.53 -15.89 800.00 0.157 0.181 -0.971
mu- 0.578 48.12 -240.97 800.00 -0.200 0.363 -0.910
mu- 2.827 -78.91 -25.19 800.00 0.024 -0.182 -0.983
   @endverbatim

   The columns are 
   - particle id:  This should be a name recognized by Gleam.
   - kinetic energy: The units, MeV or GeV, are set in the xml entry 
                     definition.
   - x, y, z: intercept point for the incident particle trajectory in 
              instrument coordinates; units are mm.
   - cos_x, cos_y, cos_z: direction cosines of the particle trajectory;
                          these need to be normalized properly.

   When the last line of particle properties is reached in the input
   file, those properties are used for all subsequent events.

   Here is an example xml definition:

   @verbatim
   <source name="file_source">
      <spectrum escale="GeV">
         <SpectrumClass name="FileSource" 
          params="input_file=$(FLUXROOT)/sources/test_FileSource.dat,rate=1,backoff_distance=0"/>
         <custom_dir/>
         <custom_pt/>
      </spectrum>
   </source>
   @endverbatim

   Note that the <tt>params</tt> string uses an interface implemented
   by Johann so that the variable names and their values are
   associated by an "=" sign. The incident particle rate is in units
   of Hz, and the "backoff_distance" is a negative offset from the
   intercept points along the particle trajectory.

   The <tt>custom_dir</tt> tag indicates that the LaunchDirection
   object provided by the FileSource object will be used instead of
   one of the those available from within the FluxSource class, e.g.,
   <tt>direction</tt>, <tt>solid_angle</tt>, <tt>celestial_dir</tt>,
   etc..  Likewise, the <tt>custom_pt</tt> tag means that the
   FileSource object's LaunchPoint object will be used instead of
   <tt>launch_point</tt> or <tt>patch</tt>.

   Either of the <tt>custom</tt> tags maybe replaced (or omitted in
   the case of <tt>custom_pt</tt>) by other valid tags, and the
   implementations indicated by those tags will be used.  For example,
   the xml entry

   @verbatim
   <source name="file_source">
      <spectrum escale="GeV">
         <SpectrumClass name="FileSource" 
          params="input_file=$(FLUXROOT)/sources/test_FileSource.dat,rate=1,backoff_distance=0"/>
         <celestial_dir ra="83.57" dec="22.01"/>
         <custom_pt/>
      </spectrum>
   </source>
   @endverbatim
   
   will result in incident particles from the direction of the Crab,
   but the specific energies and intercept points given in the input
   file will be used for each particle.
*/
