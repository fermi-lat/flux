/** @file SourceDirection.cxx
@brief SourceDirection implementation

$Header$

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
     {"zenith",  "equatorial","galactic", "Sun", "Moon"};
    if( frame == frame_names[n++] ) m_frame=ZENITH;
    if( frame == frame_names[n++] ) m_frame=EQUATORIAL;
    if( frame == frame_names[n++] ) m_frame=GALACTIC;
    if( frame == frame_names[n++] ) m_frame=SUN;
    if( frame == frame_names[n++] ) m_frame=MOON;
    if( m_frame==INVALID ){ 
        throw std::invalid_argument("flux/SourceDirection: frame name"+frame+" not recognized");
    }

}

void SourceDirection::execute(double ke, double time){
    using astro::GPS;
    using astro::SkyDir;
    using astro::SolarSystem;
    static SkyDir northPole(0,90);

    GPS* gps = GPS::instance();
    // get the direction information from the ISpectrum object
    std::pair<float,float> direction = m_spectrum->dir(ke);
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

                //here, we have a direction in the zenith direction, so we need the 
                //transformation from zenith to GLAST.
                CLHEP::HepRotation zenToGlast = gps->transformToGlast(time,GPS::ZENITH);

                setDir(zenToGlast*(-unrotated));
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

                //and do the transform, finally reversing the direction to correspond to the incoming particle
                setDir( - (celtoglast * unrotated()) );
                break;
            }
 
        case SUN:
            {
                // get the direction to the Sun or Moon
  
                // direction to center of the object, polar angle
                Hep3Vector dir(solarsystemDir(time));
                double theta(dir.theta());

                // direction with respect to object
                Hep3Vector axis(northPole().cross(dir));
                double costh = first
                    ,  sinth = sqrt(1.-costh*costh)
                    ,  phi   = second;
                CLHEP::Hep3Vector unrotated(cos(phi)*sinth, sin(phi)*sinth, costh);

                HepRotation R(axis, theta);
                setDir(- (R* unrotated));
                break;
            }
        case MOON:
            {
                SolarSystem luna(astro::SolarSystem::MOON);
                double jd(astro::JulianDate::missionStart()+time);

                setDir(-luna.direction(jd).dir());    
                break;
            }
    }

}

Hep3Vector solarsytemDir( double time)
{

    double jd(astro::JulianDate::missionStart()+time);

    if( m_frame==SUN) {
        SolarSystem sol(astro::SolarSystem::SUN);
        return Hep3Vector dir(sol.direction(jd)());
    }
    if (m_frame==MOON) {
        SolarSystem luna(astro::SolarSystem::MOON);

        return Hep3Vector dir(luna.direction(jd)());


    }
}


//! solid angle
double SourceDirection::solidAngle()const {
    return m_spectrum->solidAngle();
}

