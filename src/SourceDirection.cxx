/** @file SourceDirection.cxx
@brief SourceDirection implementation

$Header: /nfs/slac/g/glast/ground/cvs/flux/src/SourceDirection.cxx,v 1.1 2006/12/03 03:36:08 burnett Exp $

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
        case MOON:
            {
                solarSystemDir(first, second, time);
                break;
            }
    }

}

void SourceDirection::solarSystemDir( double costh, double phi, double time)
{
    using astro::SolarSystem;
    using CLHEP::Hep3Vector;
    static Hep3Vector northPole(0,0,1);

    double  sinth( sqrt(1.-costh*costh));
    Hep3Vector unrotated(cos(phi)*sinth, sin(phi)*sinth, costh);
    double jd(astro::JulianDate::missionStart()+time);

    Hep3Vector dir;
    if( m_frame==SUN) {
        SolarSystem sol(astro::SolarSystem::SUN);
        dir = Hep3Vector(sol.direction(jd)());
    }
    if (m_frame==MOON) {
        SolarSystem luna(astro::SolarSystem::MOON);
        dir = Hep3Vector(luna.direction(jd)());
    }
    Hep3Vector axis( northPole.cross(dir));
    HepRotation R(axis, dir.theta());

    setDir(- (R* unrotated));

}


//! solid angle
double SourceDirection::solidAngle()const {
    return m_spectrum->solidAngle();
}

