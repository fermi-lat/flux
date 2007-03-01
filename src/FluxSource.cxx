/** @file FluxSource.cxx
@brief Implementation of FluxSource

$Header: /nfs/slac/g/glast/ground/cvs/flux/src/FluxSource.cxx,v 1.42 2006/12/22 20:31:59 burnett Exp $

*/
#include "astro/SkyDir.h"
#include "astro/GPS.h"

#include "flux/FluxSource.h"
#include "flux/LaunchDirection.h"
#include "flux/LaunchPoint.h"
#include "flux/SpectrumFactoryTable.h"
#include "flux/SimpleSpectrum.h"
#include "flux/FluxException.h" // for FATAL_MACRO

#include "SourceDirection.h"

#include <xercesc/dom/DOMElement.hpp>
#include "xmlBase/Dom.h"

#include "CLHEP/Random/RandFlat.h"


#include <algorithm>
#include <sstream>
#include <stdexcept>


//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
/** @class RandomPoint
@brief nested launch strategy derived class
This is the standard strategy, which takes a direction and creates a point in
a disk centered at the origin with area 6 m^2 (or so)
*/
class FluxSource::RandomPoint : public LaunchPoint { 
public:
    RandomPoint(double radius, double backoff)
        :m_radius(radius), m_backoff(backoff)
    { 

    }

    virtual void execute(const CLHEP::Hep3Vector& dir){
        CLHEP::HepRotation r_pln;

        //create rotation to take x-y plane to be perpendicular to incoming direction
        double ly = dir.y(), lx = dir.x();
        if( fabs( lx) +fabs(ly) >1e-8) {  // leave as identity 
            r_pln.rotate(acos(dir.z()),  CLHEP::Hep3Vector(-ly, lx, 0.));
        }

        // pick a random position on the planar section of a sphere through 
        // its midpoint
        double 
            azimuth = CLHEP::RandFlat::shoot( 2*M_PI ),
            rad = m_radius*(sqrt(CLHEP::RandFlat::shoot()));

        // create two vectors to describe the particle launch: one to describe
        // the point in the plane perpendicular to the launch direction (within
        // the cross-section of the sphere containing the instrument) and a 
        // second to describe the distance along the normal between the launch 
        // point and that plane.
#if 1 // standard
        CLHEP::Hep3Vector posLaunch(rad*cos(azimuth), rad*sin(azimuth), 0.);

        // define actual launch point
        setPoint( r_pln*posLaunch - m_backoff*dir);
#else // patch for Atwood calorimeter study
        HepPoint3D posLaunch(RandFlat::shoot(180.), 180., -56.);
	setPoint( posLaunch - m_backoff*dir);
#endif
    }

    /// return info, 
    std::string title()const{
        std::stringstream t;
        t << "radius("<< m_radius << ")";
        return t.str();
    }


private:
    double m_radius;
    double m_backoff;

}; 
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
/** @class FixedPoint
@brief nested launch strategy derived class

This strategy uses a fixed launch point for a pencil beam. If the radius is nonzero,
the beam will be spread out uniformly on a disk perpendicular to the incoming direction
*/
class FluxSource::FixedPoint : public LaunchPoint{ 
public:
    FixedPoint( const CLHEP::Hep3Vector& pt, double radius)
        :  LaunchPoint(pt)
        ,  m_disk_radius(radius)
        ,  m_base_point(pt)
    {}

    virtual void execute(const CLHEP::Hep3Vector& dir){
        if(m_disk_radius==0) return; // just use 

        CLHEP::HepRotation r_pln;

        double ly = dir.y(), lx = dir.x();
        if( lx !=0 || ly !=0 ) { 
            r_pln.rotate(acos(dir.z()), CLHEP::Hep3Vector(-ly, lx, 0.));
        }
        double 
            azimuth = CLHEP::RandFlat::shoot( 2*M_PI ),
            rad = m_disk_radius*(sqrt(CLHEP::RandFlat::shoot()));
        CLHEP::Hep3Vector posLaunch(rad*cos(azimuth), rad*sin(azimuth), 0.);

        setPoint(r_pln*posLaunch + m_base_point);

    };
    virtual std::string title() const {
        if( m_disk_radius==0) return LaunchPoint::title();
        std::stringstream t;
        t << ", radius(" << m_disk_radius << ")";
        return LaunchPoint::title() + t.str();
    }
private:
    double m_disk_radius;
    CLHEP::Hep3Vector m_base_point;
};  
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
/** @class Patch
@brief nested launch strategy derived class
Gets a point randomly from a box
*/
class FluxSource::Patch : public LaunchPoint { 
public:
    Patch( double xmin, double xmax, double ymin, double ymax, double zmin, double zmax)
        :m_xmin(xmin), m_dx(xmax-xmin), 
        m_ymin(ymin), m_dy(ymax-ymin), 
        m_zmin(zmin), m_dz(zmax-zmin)
    {
    }

    virtual void execute(const CLHEP::Hep3Vector& ){
        setPoint(CLHEP::Hep3Vector( 
            m_xmin + m_dx*CLHEP::RandFlat::shoot(),
            m_ymin + m_dy*CLHEP::RandFlat::shoot(),
            m_zmin + m_dz*CLHEP::RandFlat::shoot()) );
    }
    virtual std::string title() const {
        std::stringstream t;
        t << "patch(" 
            << m_xmin << "," << m_xmin+m_dx << ","
            << m_ymin << "," << m_ymin+m_dy << ","
            << m_zmin << "," << m_zmin+m_dz << ")" ;
        return t.str();
    }

private:
    double m_xmin, m_dx, m_ymin, m_dy, m_zmin, m_dz;    
}; 
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
/** @class RandomDirection
@brief nested launch strategy derived class
Assigns a random direction from a range of cos theta, optionally rotated
*/
class FluxSource::RandomDirection : public LaunchDirection{ 
public:
    /** ctor:
    @param minc  minimum value of cos(theta)
    @param maxc  maximum value of cos(theta)
    @param theta [0] X rotation angle (radians)
    @param phi   [0] Z rotation angle (radians)
    */
    RandomDirection(double minc, double maxc, double theta=0, double phi=0)
        : m_theta(theta)
        , m_phi(phi)
    {
        using std::min;
        using std::max;

        // require _maxCos > _minCos for solid angle calculation
        m_minCos   = min(minc,maxc);
        m_maxCos   = max(minc,maxc);
        if(m_minCos==m_maxCos) {
            if(m_minCos!=-1) m_minCos-=0.001; else m_maxCos +=0.001;
        }
        m_minPhi = 0; m_maxPhi=2*M_PI;

    }

   virtual void execute(double /*ke*/, double time){
        double  costh = -CLHEP::RandFlat::shoot(m_minCos, m_maxCos),
            sinth = sqrt(1.-costh*costh),
            phi = CLHEP::RandFlat::shoot(m_minPhi, m_maxPhi);

        //here, the direction is with respect to the zenith frame,
        //so we need the transformation from the zenith to GLAST.
        CLHEP::HepRotation zenToGlast=astro::GPS::instance()->transformToGlast(time,astro::GPS::ZENITH);

        CLHEP::Hep3Vector dir(cos(phi)*sinth, sin(phi)*sinth, costh);

        // extra rotation in case not zenith pointing (beware, might be
        // confusing)
        // keep x-axis perpendicular to zenith direction
        if (m_theta != 0.0) dir.rotateX(m_theta).rotateZ(m_phi);
        setDir(zenToGlast*dir);

    }
    //! solid angle
    virtual double solidAngle()const {
        return 2*M_PI*(m_maxCos-m_minCos);
    }


    virtual std::string title()const {
        std::stringstream t;
        t << "range(" << m_minCos << ',' << m_maxCos << ") ";
        if( m_theta != 0){
            t << ", angle(" << m_theta*180/M_PI << ',' << m_phi*180/M_PI << ") ";
        }
        return t.str();
    }


private:
    double m_minCos, m_maxCos;
    double m_minPhi, m_maxPhi;
    double m_theta, m_phi;

}; 

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//                   FluxSource constructor
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
FluxSource::FluxSource(const XERCES_CPP_NAMESPACE_QUALIFIER DOMElement* xelem )
: EventSource ()
, m_spectrum(0)
, m_occultable(true)
, m_zenithCosTheta(1.0) //won't be occulted by default
, m_launch_dir_owner(true)
, m_launch_pt_owner(true)

{
    using XERCES_CPP_NAMESPACE_QUALIFIER DOMElement;
    static double d2r = M_PI/180.;
    setName(xmlBase::Dom::getAttribute(xelem, "name").c_str());


    ISpectrum*   s = 0;
    std::string class_name;
    std::string source_params; 
    // this is a default flux, from the flux="123" in the source element
    setFlux(atof (xmlBase::Dom::getAttribute(xelem, "flux").c_str()));

    DOMElement*   spec = xmlBase::Dom::findFirstChildByName(xelem, "spectrum");

    if (spec == 0) {

        // source has no imbedded spectrum element: expect a name
        class_name = xmlBase::Dom::getAttribute(xelem, "name");
    } else {
        // process spectrum element
        DOMElement* specType = xmlBase::Dom::getFirstChildElement(spec);

        std::string typeTagName = xmlBase::Dom::getTagName(specType);
        std::string particle_name = xmlBase::Dom::getAttribute(spec, "particle_name");
        std::string spectrum_energyscale = xmlBase::Dom::getAttribute(spec, "escale");

        if(spectrum_energyscale == "GeV"){ m_energyscale=GeV;
        }else if(spectrum_energyscale == "MeV"){ m_energyscale=MeV;
        }else{
            std::cout << "bad energy scale declaration on spectrum:"
                << spectrum_energyscale << " , exiting.";
            return;} //this line "just in case"

        if (typeTagName=="particle") s = new SimpleSpectrum(specType,  m_energyscale==GeV );
        else if (typeTagName=="SpectrumClass") {
            // attribute "name" is the class name
            class_name = xmlBase::Dom::getAttribute(specType, "name");
            source_params= xmlBase::Dom::getAttribute(specType, "params");
        }
        else {
            // no, the tag itself
            class_name = typeTagName;//.transcode();
        }

        //if s is still 0, we need to create the internal spectrum object.
        if( s==0) {
            //		std::vector<float> paramvec; parseParamList(source_params, paramvec);
            s = SpectrumFactoryTable::instance()->instantiate(class_name, source_params);
            if(s==0){

                std::cerr << "List of known Spectrum classes:\n" ;
                std::list<std::string>list= SpectrumFactoryTable::instance()->spectrumList();
                for( std::list<std::string>::iterator i = list.begin(); i!=list.end(); ++i)
                    std::cerr << "\t" << *i << std::endl;
                FATAL_MACRO("Unknown Spectrum: "<< class_name);
                return;
            }
	    std:: string flux = xmlBase::Dom::getAttribute(spec, "flux");
	    if(flux!="")
	      { s->setFlux(atof(flux.c_str())); }
	    s->setInGeV(spectrum_energyscale == "GeV");

	    if( !particle_name.empty() ) s->setParticleName(particle_name);
        }
        m_spectrum =s;
	

        // second child element is angle
        DOMElement* angles = xmlBase::Dom::getSiblingElement(specType);
        std::string anglesTag = xmlBase::Dom::getTagName(angles);
        if (anglesTag == "solid_angle") 
        {
            m_occultable=false;
            m_launch_dir = new RandomDirection(
                xmlBase::Dom::getDoubleAttribute(angles, "mincos"),
                xmlBase::Dom::getDoubleAttribute(angles, "maxcos"),
                xmlBase::Dom::getDoubleAttribute(angles, "theta") * d2r,
                xmlBase::Dom::getDoubleAttribute(angles, "phi")*d2r);

        }
        else if (anglesTag == "direction") 
        {
            //m_occultable=false;
            std::string frame = xmlBase::Dom::getAttribute(angles, "frame");
            m_occultable=(frame=="zenith");
            m_launch_dir = new LaunchDirection(
                xmlBase::Dom::getDoubleAttribute(angles, "theta") * d2r,
                xmlBase::Dom::getDoubleAttribute(angles, "phi")*d2r,
                frame);
        }
        else if (anglesTag == "use_spectrum")
        {
            std::string frame = 
                xmlBase::Dom::getAttribute(angles, "frame");
            m_occultable=(frame!="zenith");
            m_launch_dir = new SourceDirection(m_spectrum, frame); 
        }
        else if (anglesTag == "galactic_dir")
        {
            m_occultable=true;
            m_launch_dir = new LaunchDirection(
                astro::SkyDir(
                xmlBase::Dom::getDoubleAttribute(angles, "l") ,
                xmlBase::Dom::getDoubleAttribute(angles, "b") ,
                astro::SkyDir::GALACTIC
                ),
               xmlBase::Dom::getDoubleAttribute(angles, "radius") 

                );
        }
        else if (anglesTag == "celestial_dir")

        {
            m_occultable=true;
            m_launch_dir = new LaunchDirection(
                astro::SkyDir(
                xmlBase::Dom::getDoubleAttribute(angles, "ra")  ,
                xmlBase::Dom::getDoubleAttribute(angles, "dec")

                ),
               xmlBase::Dom::getDoubleAttribute(angles, "radius") 
                );

        }
        else if (anglesTag == "galactic_spread")
        {
            m_occultable=true;
            FATAL_MACRO("not implemented");
        }
        else if (anglesTag == "custom_dir") {
// These sources are not intended for modeling astrophysical objects
// and so cannot be occulted.
           m_occultable = false;     
           m_launch_dir = m_spectrum->launchDirection();
           if (m_launch_dir == 0) {
              std::ostringstream what;
              what << "FluxSource: cannot use a 'custom_dir' tag with a "
                   << class_name << " source.";
              throw std::runtime_error(what.str());
           }
           m_launch_dir_owner = false;
        }
        else {
            FATAL_MACRO("Unknown angle specification in Flux::Flux \""
                << anglesTag << "\"" );
        }

        // third child element is optional launch spec
        DOMElement* launch = xmlBase::Dom::getSiblingElement(angles);

        if(launch != 0) {
            std::string launchTag = xmlBase::Dom::getTagName(launch);

            if (launchTag == "launch_point")
            {
                m_launch_pt = new FixedPoint(CLHEP::Hep3Vector(
                    xmlBase::Dom::getDoubleAttribute(launch, "x"),
                    xmlBase::Dom::getDoubleAttribute(launch, "y"),
                    xmlBase::Dom::getDoubleAttribute(launch, "z")),
                    xmlBase::Dom::getDoubleAttribute(launch, "beam_radius") );
            }
            else if (launchTag == "patch")
            {
                m_launch_pt = new Patch(
                    xmlBase::Dom::getDoubleAttribute(launch, "xmax"),
                    xmlBase::Dom::getDoubleAttribute(launch, "xmin"),
                    xmlBase::Dom::getDoubleAttribute(launch, "ymax"),
                    xmlBase::Dom::getDoubleAttribute(launch, "ymin"),
                    xmlBase::Dom::getDoubleAttribute(launch, "zmax"),
                    xmlBase::Dom::getDoubleAttribute(launch, "zmin") );
            } else if (launchTag == "custom_pt") {
               m_launch_pt = m_spectrum->launchPoint();
               if (m_launch_pt == 0) {
                  std::ostringstream what;
                  what << "FluxSource: cannot use a 'custom_pt' tag with a "
                       << class_name << " source.";
                  throw std::runtime_error(what.str());
               }
               m_launch_pt_owner = false;
            } else {
               FATAL_MACRO("Unknown launch specification in Flux::Flux \""
                           << launchTag << "\"" );
            }
        } else {
            // default: the target sphere.
            double radius = sqrt(totalArea() / M_PI ) * 1000;   // radius in mm
            m_launch_pt = new RandomPoint(radius, EventSource::s_backoff);
        }
    }
}


FluxSource::~FluxSource()
{
    delete m_spectrum;
    if (m_launch_pt_owner) delete m_launch_pt;
    if (m_launch_dir_owner) delete m_launch_dir;
}

void FluxSource::spectrum(ISpectrum* s, double emax)
{
    if (emax > 0) {
        //m_rmax =  s->fraction(emax);
        std::cerr << "exercising obsolete function fraction" << std::endl;
    }
    m_spectrum = s;
    //const char* name = s->particleName();
}

EventSource* FluxSource::event(double time)
{
    // Purpose and Method: generate a new incoming particle
    // Inputs  - current time
    // Outputs - pointer to the "current" fluxSource object,or zero if it has "turned off"
    if( !enabled()){
        throw std::runtime_error("FluxSource::event called when disabled");
    }
    using astro::GPS;
    m_interval = calculateInterval(time);
    if( time+m_interval < GPS::instance()->endTime()){
        // do this only if in valid interval: assume will never get used otherwise
        computeLaunch(time + m_interval);
    }else{
        // flag to end use of this source
       disable();
    }
    //now set the actual interval to be what FluxMgr will get, unless beyond the endtime
    EventSource::setTime(time + m_interval);
    return this;
}

double FluxSource::calculateInterval (double time)
{   
    double interval=m_spectrum->interval(time);
    if( interval>0 ){
        // the spectum computed an interval: use it
        return interval;
    }

    // otherwise do a Poison from the the flux, solid angle, and area factor
    return -log(1.-CLHEP::RandFlat::shoot(1.))/rate(time);
}

double FluxSource::flux(double time) const
{
    if(!enabled()){ return 0;}
    if(m_spectrum->flux(time)){ return m_spectrum->flux(time);}
    else{return EventSource::flux(time);}
}

double FluxSource::rate(double time) const
{
    //TODO: area is only relevant for RandomPoint, and flux per unit area
    return m_launch_dir->solidAngle() * flux(time) * totalArea();
}

void FluxSource::computeLaunch (double time)
{
    // get the KE from the spectrum object
    m_energy = spectrum()->energy( time );

    // convert to MeV if necessary
    if(m_energyscale==GeV){
        m_energy *= 1000.;
    }

    // set launch direction , position (perhaps depending on direction)
    m_launch_dir->execute(m_energy, time);
    m_launchDir  = (*m_launch_dir)();

    //get the off-zenith angle cosine, for occultation purposes:
    m_zenithCosTheta = m_launch_dir->zenithCosine();

    //  rotate by Glast orientation transformation
    if( s_applyAlign) m_correctedDir = s_alignMatrix * m_launchDir;
    else     m_correctedDir = m_launchDir;

    // now set the launch point, which may depend on the direction

    m_launch_pt->execute(m_correctedDir);

    m_launchPoint = (*m_launch_pt)();

}

std::string FluxSource::fullTitle () const
{
    return title();
}

std::string FluxSource::displayTitle () const
{
    std::stringstream s;
    s << EventSource::displayTitle() << '(' << m_spectrum->title() ;
    s << ')' << '\0';
    return s.str();
}

int FluxSource::eventNumber()const
{
    return 0;
}


std::string FluxSource::title () const
{
    return m_spectrum->title() + ", "
        +  m_launch_pt->title() +", "
        +  m_launch_dir->title();
}

std::string FluxSource::particleName()
{
    return spectrum()->particleName();
}


bool FluxSource::occulted(){
    //Purpose:  to determine whether or not the current incoming particle will be blocked by the earth.
    //Output:  "yes" or "no"
    //REMEMBER:  the earth is directly below the satellite, so, to determine occlusion,
    // we must assume the frame to be checked against is zenith-pointing, and hence, we want 
    //the direction of the particle relative to the GLAST zenith direction, calculated in the 
    //LaunchDirection classes.

    //this should probably be open-ended, not wired in!
    double minCosTheta= -0.4;

    return (m_occultable) && (m_zenithCosTheta < minCosTheta);

}

astro::SkyDir FluxSource::skyDirection()const
{
    return astro::SkyDir(m_launch_dir->dir());
}

void FluxSource::disable()
{
    m_enabled=false;
    delete m_spectrum;
    m_spectrum=0;
}

