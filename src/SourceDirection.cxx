/** @file SourceDirection.cxx
@brief SourceDirection implementation

$Header: /nfs/slac/g/glast/ground/cvs/ScienceTools-scons/flux/src/SourceDirection.cxx,v 1.14 2011/05/20 16:14:56 heather Exp $

*/

#include "SourceDirection.h"
#include "flux/ISpectrum.h"

#include "astro/SolarSystem.h"
#include "astro/JulianDate.h"

#include <vector>
#include <stdexcept>

SourceDirection::SourceDirection(ISpectrum* spectrum, std::string frame )
: m_spectrum(spectrum)
, m_frameName(frame)
, m_zenithCos(1.0)

{
    m_frame = INVALID; int n(0);
    static const char* frame_names[]=
     {"zenith",  "equatorial","galactic", "galaxy", "Sun", "Moon", "Jupiter", "Saturn"};
    if( frame == frame_names[n++] ) m_frame=ZENITH;
    if( frame == frame_names[n++] ) m_frame=EQUATORIAL;
    if( frame == frame_names[n++] ) m_frame=GALACTIC;
    if( frame == frame_names[n++] ) m_frame=GALACTIC;
    if( frame == frame_names[n++] ) m_frame=SUN;
    if( frame == frame_names[n++] ) m_frame=MOON;
    if( frame == frame_names[n++] ) m_frame=JUPITER;
    if( frame == frame_names[n++] ) m_frame=SATURN;
    if( m_frame==INVALID ){ 
        throw std::invalid_argument("flux/SourceDirection: frame name"+frame+" not recognized");
    }
    

}

void SourceDirection::execute(double ke, double time){
    using astro::GPS;
    using astro::SkyDir;
    using astro::SolarSystem;

    GPS* gps = GPS::instance();
    // get the direction information from the ISpectrum object
    std::pair<double,double> direction = m_spectrum->dir(ke);
    double first(direction.first), second(direction.second);

    switch (m_frame) {
        case ZENITH:
            {
                // special option that gets direction from the spectrum object
                // note extra - sign since direction corresponds to *from*, not *to*

                double  costh = first,
                    sinth = sqrt(1.-costh*costh),
                    phi = second;
                CLHEP::Hep3Vector unrotated(cos(phi)*sinth, sin(phi)*sinth, costh);

                if (::getenv("ZENITH_FRAME_FIX")) {
                   // Use right-handed coordinate system with z-axis
                   // along nadir and x-axis to North
                   unrotated = CLHEP::Hep3Vector(cos(phi)*sinth, sin(phi)*sinth, -costh);
                }
                
                //here, we have a direction in the zenith direction, so we need the 
                //transformation from zenith to GLAST.
                CLHEP::HepRotation zenToGlast = gps->transformToGlast(time,GPS::ZENITH);

                // use new transformation
                m_lat_dir = - gps->LATdirection(GPS::ZENITH, unrotated) ;

                break;
            }
        case EQUATORIAL:
        case GALACTIC:
            {
                // interpret direction as (l,b) or (ra,dec) for a  celestial source
                //then set up this direction, either in galactic or equatorial coordinates:    
                astro::SkyDir unrotated(first, second,
                    m_frame==GALACTIC? SkyDir::GALACTIC : SkyDir::EQUATORIAL);

                //get the zenith cosine:
                astro::SkyDir zenDir(gps->zenithDir());
                m_zenithCos = unrotated()*zenDir();
                //get the transformation matrix..
                //here, we have a SkyDir, so we need the transformation from a SkyDir to GLAST.
                CLHEP::HepRotation celtoglast = gps->transformToGlast(time, GPS::CELESTIAL);

                // transform from sky direction to LAT direction, possibly with aberration and alignment
                m_lat_dir = gps->LATdirection(GPS::CELESTIAL, unrotated()) ;

                break;
            }
 
        case SUN:
        case MOON:
        case JUPITER:
        case SATURN:
            {
                solarSystemDir(first, second, time);
                break;
            }
    }

}

void SourceDirection::solarSystemDir( double ra, double dec, double time)
{
    // expect displacement with respect to the object's direction

    using astro::GPS;    
    using astro::SolarSystem;
    using astro::SkyDir;
    using astro::JulianDate;
    using CLHEP::Hep3Vector;

    GPS* gps = GPS::instance();
    static CLHEP::Hep3Vector xhat(1,0,0);

    JulianDate jd(JulianDate::missionStart()+time/JulianDate::secondsPerDay);

    Hep3Vector cdir;
    if( m_frame==SUN) {
        SolarSystem sol(astro::SolarSystem::SUN);
        cdir = CLHEP::Hep3Vector(sol.direction(jd)());
    }
    if( m_frame==JUPITER) {
        cdir = CLHEP::Hep3Vector(SolarSystem(SolarSystem::JUPITER).direction(jd)());
    }
    if( m_frame==SATURN) {
        cdir = CLHEP::Hep3Vector(SolarSystem(SolarSystem::SATURN).direction(jd)());
    }
    if (m_frame==MOON) {
        SolarSystem luna(astro::SolarSystem::MOON);
        cdir = CLHEP::Hep3Vector(luna.direction(jd, gps->position())());
    }
    // create rotation that takes (0,0) to (ra,dec)
    CLHEP::Hep3Vector input ( SkyDir(ra,dec)() );
    double  theta ( acos(input.x()) )
        ,   phi ( atan2(input.y(), input.z()) );

    // first an axis perpendicular to the given direction
    CLHEP::Hep3Vector axis( cdir.cross(xhat).unit() );
    if(cdir.isNear(xhat)) axis = CLHEP::Hep3Vector(0,1,0);

    // rotate that axis by phi about the solar object direction
    CLHEP::Hep3Vector axis_prime ( CLHEP::HepRotation(cdir, phi) * axis );

    // and a rotation about the new axis by theta
    CLHEP::Hep3Vector rdir ( CLHEP::HepRotation(axis_prime,theta) * cdir );

#if 0 // stuff for debugging checks
    double open( acos(rdir*cdir)*180/M_PI); // should be the total rotation
    SkyDir r(rdir); double ra_r(r.ra()), dec_r(r.dec());
    SkyDir c(cdir); double ra_c(c.ra()), dec_c(c.dec());

#endif

    // make transformation, to a direction into the LAT, possibly 
    // corrected for misalignment and aberration
    m_lat_dir = gps->LATdirection( GPS::CELESTIAL, rdir , time);
}


//! solid angle
double SourceDirection::solidAngle()const {
    return m_spectrum->solidAngle();
}

