/** @file FluxSource.cxx
@brief Implementation of FluxSource

$Header: /nfs/slac/g/glast/ground/cvs/flux/src/FluxSource.cxx,v 1.13 2004/01/15 12:27:03 burnett Exp $

*/
#include "flux/FluxSource.h"

#include <xercesc/dom/DOM_Element.hpp>
#include <xercesc/dom/DOM_NodeList.hpp>
#include "xml/Dom.h"
#include "CLHEP/Random/RandFlat.h"

#include "astro/SkyDir.h"

#include "flux/SpectrumFactoryTable.h"
#include "flux/SimpleSpectrum.h"

#include "flux/FluxException.h" // for FATAL_MACRO
#include "flux/GPS.h"

#include <algorithm>
#include <sstream>
namespace {
    // this is the (wired-in) distance to back off from the target sphere.
    const double backoff_distance=1500.;
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
/** @class LaunchPoint
@brief nested launch strategy base class for point determination

The virtual base class manages the point itself
*/
class FluxSource::LaunchPoint  { 
public:
    LaunchPoint(){}
    LaunchPoint(const HepPoint3D& pt):m_pt(pt){}
    virtual ~LaunchPoint(){}

    /// access to direction, perhaps set by the execute()
    virtual const HepPoint3D& point()const {return m_pt;}
    const HepPoint3D& operator()()const{return point();}

    /// execute the strategy, perhaps depending on direction
    virtual void execute(const HepVector3D& ){};

    /// set the point
    void setPoint(const HepPoint3D& pt){ m_pt = pt;}

    /// return info, default if not overriden
    virtual std::string title()const{
        std::stringstream t;
        t << "point" << m_pt;
        return t.str();
    }

private:
    HepPoint3D m_pt;
};
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
/** @class RandomPoint
@brief nested launch strategy derived class
This is the standard strategy, which takes a direction and creates a point in
a disk centered at the origin with area 6 m^2 (or so)
*/
class FluxSource::RandomPoint : public LaunchPoint{ 
public:
    RandomPoint(double radius, double backoff)
        :m_radius(radius), m_backoff(backoff)
    { 

    }

    virtual void execute(const HepVector3D& dir){
        HepRotation r_pln;

        //create rotation to take x-y plane to be perpendicular to incoming direction
        double ly = dir.y(), lx = dir.x();
        if( fabs( lx) +fabs(ly) >1e-8) {  // leave as identity 
            r_pln.rotate(acos(dir.z()),  HepVector3D(-ly, lx, 0.));
        }

        // pick a random position on the planar section of a sphere through 
        // its midpoint
        double 
            azimuth = RandFlat::shoot( 2*M_PI ),
            rad = m_radius*(sqrt(RandFlat::shoot()));

        // create two vectors to describe the particle launch: one to describe
        // the point in the plane perpendicular to the launch direction (within
        // the cross-section of the sphere containing the instrument) and a 
        // second to describe the distance along the normal between the launch 
        // point and that plane.
        HepPoint3D posLaunch(rad*cos(azimuth), rad*sin(azimuth), 0.);

        // define actual launch point
        setPoint( r_pln*posLaunch - m_backoff*dir);
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
    FixedPoint( const HepPoint3D& pt, double radius)
        :  LaunchPoint(pt)
        ,  m_disk_radius(radius)
        ,  m_base_point(pt)
    {}

    virtual void execute(const HepVector3D& dir){
        if(m_disk_radius==0) return; // just use 

        HepRotation r_pln;

        double ly = dir.y(), lx = dir.x();
        if( lx !=0 || ly !=0 ) { 
            r_pln.rotate(acos(dir.z()), HepVector3D(-ly, lx, 0.));
        }
        double 
            azimuth = RandFlat::shoot( 2*M_PI ),
            rad = m_disk_radius*(sqrt(RandFlat::shoot()));
        HepPoint3D posLaunch(rad*cos(azimuth), rad*sin(azimuth), 0.);

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
    HepPoint3D m_base_point;
};  
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
/** @class Patch
@brief nested launch strategy derived class
Gets a point randomly from a box
*/
class FluxSource::Patch : public FluxSource::LaunchPoint{ 
public:
    Patch( double xmin, double xmax, double ymin, double ymax, double zmin, double zmax)
        :m_xmin(xmin), m_dx(xmax-xmin), 
        m_ymin(ymin), m_dy(ymax-ymin), 
        m_zmin(zmin), m_dz(zmax-zmin)
    {
    }

    virtual void execute(const HepVector3D& ){
        setPoint(HepPoint3D( 
            m_xmin + m_dx*RandFlat::shoot(),
            m_ymin + m_dy*RandFlat::shoot(),
            m_zmin + m_dz*RandFlat::shoot()) );
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
/** @class LaunchDirection
@brief nested launch strategy base class
*/
class FluxSource::LaunchDirection  {
public:
    LaunchDirection():m_skydir(false),m_radius(0){}

    virtual ~LaunchDirection(){}

    LaunchDirection(double theta, double phi, double radius=0)
        :m_skydir(false)
        , m_radius(radius*M_PI/180)
    {
        HepVector3D dir(sin(theta)*cos(phi),sin(theta)*sin(phi),cos(theta));
        setDir(-dir); // minus due to z axis pointing UP!
    }
    LaunchDirection(astro::SkyDir sky, double radius=0)
        :m_skydir(true)
        , m_radius(radius*M_PI/180)
    {
        m_dir = -sky.dir();
    }
    /** @brief choose a direction
    @param KE kinetic energy
    @param time mission time
    */
    virtual void execute(double KE, double time){
        if(m_skydir){
            //here, we have a SkyDir, so we need the transformation from a SkyDir to GLAST.
            m_rottoglast = GPS::instance()->transformToGlast(time,GPS::CELESTIAL);//->transformCelToGlast(time);
        }else{
            //otherwise, the direction is in the zenith system, and the rotation to GLAST is needed:
            m_rottoglast = GPS::instance()->transformToGlast(time,GPS::ZENITH);
        }
    }

    const HepVector3D& operator()()const {return dir();}

    virtual const HepVector3D& dir()const {
        static HepVector3D rdir;
        rdir = m_rottoglast * m_dir;
        if( m_radius>0 ) {
            // spread uniformly about a disk
            // rotate about perpendicular then about the original 
            HepVector3D t(rdir);
            t.rotate( m_radius*(sqrt(RandFlat::shoot())),  rdir.orthogonal()),  // rotate about the orthogonal
            t.rotate( RandFlat::shoot( 2*M_PI ), rdir); // rotate about the original direction
            rdir = t; //replace 
        }
        return rdir;
    }

    void setDir(const HepVector3D& dir){m_dir=dir;}


    //! solid angle: default of 1. for a point source
    virtual double solidAngle()const {
        return 1.;
    }

    /// return info, default if not overriden
    virtual std::string title()const{
        std::stringstream t;
        t << " dir" << m_dir ;
        if( m_radius>0 ) { t << " radius " << m_radius ;}
        return t.str();
    }

    /// return the cosine of the angle between the incoming direction and the earth's zenith
    virtual double zenithCosine()const{
        if(m_skydir){
            astro::SkyDir zenDir(GPS::instance()->RAZenith(),GPS::instance()->DECZenith());
            return -m_dir*zenDir();
        }
        //if the direction is local
        return 1.0;
    }

private:
    HepRotation m_rottoglast;
    HepVector3D m_dir;
    bool  m_skydir;
    HepVector3D m_t;
    double m_radius;

};
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
/** @class RandomDirection
@brief nested launch strategy derived class
Assigns a random direction from a range of cos theta, optionally rotated
*/
class FluxSource::RandomDirection : public FluxSource::LaunchDirection{ 
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

    virtual void execute(double ke, double time){
        double  costh = -RandFlat::shoot(m_minCos, m_maxCos),
            sinth = sqrt(1.-costh*costh),
            phi = RandFlat::shoot(m_minPhi, m_maxPhi);

        //here, the direction is with respect to the zenith frame,
        //so we need the transformation from the zenith to GLAST.
        HepRotation zenToGlast=GPS::instance()->transformToGlast(time,GPS::ZENITH);

        HepVector3D dir(cos(phi)*sinth, sin(phi)*sinth, costh);

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
/** @class SourceDirection
@brief nested launch strategy derived class
Gets a direction from the ISpectrum class
*/
class FluxSource::SourceDirection : public FluxSource::LaunchDirection{ 
public:
    /** Ctor:
    @param spectrum pointer to the ISpectrum object that will provide the direction
    @param galactic if true, interpret pair as l,b (in degrees); otherwise costh, phi
    */
    SourceDirection(ISpectrum* spectrum, std::string frame /*bool galactic*/)
        : m_spectrum(spectrum)
        , m_galactic(frame=="galaxy")
        , m_equatorial(frame=="equatorial")
        , m_zenithCos(1.0)
    {}

    void execute(double ke, double time){

        std::pair<float,float> direction 
            = m_spectrum->dir(ke);

        if( !(m_galactic||m_equatorial) ) {
            // special option that gets direction from the spectrum object
            // note extra - sign since direction corresponds to *from*, not *to*

            double  costh = direction.first,
                sinth = sqrt(1.-costh*costh),
                phi = direction.second;

            //here, we have a direction in the zenith direction, so we need the 
            //transformation from zenith to GLAST.
            HepRotation zenToGlast = GPS::instance()->transformToGlast(time,GPS::ZENITH);

            HepVector3D unrotated(cos(phi)*sinth, sin(phi)*sinth, costh);

            //setDir(-HepVector3D(cos(phi)*sinth, sin(phi)*sinth, costh));
            setDir(zenToGlast*(-unrotated));

        }else {
            // iterpret direction as l,b for a galactic or celestial source
            double  l = direction.first,
                b = direction.second;
            //then set up this direction, either in galactic or celestial coordinates:    
            astro::SkyDir unrotated(l,b,m_galactic? astro::SkyDir::GALACTIC : astro::SkyDir::CELESTIAL);

            //get the zenith cosine:
            astro::SkyDir zenDir(GPS::instance()->RAZenith(),GPS::instance()->DECZenith());
            m_zenithCos = unrotated()*zenDir();
            //get the transformation matrix..
            //here, we have a SkyDir, so we need the transformation from a SkyDir to GLAST.
            HepRotation celtoglast = GPS::instance()->transformToGlast(time,GPS::CELESTIAL);

            //and do the transform, finally reversing the direction to correspond to the incoming particle
            setDir( - (celtoglast * unrotated()) );
        }
    }

    //! solid angle
    virtual double solidAngle()const {
        return m_spectrum->solidAngle();
    }

    virtual std::string title()const {
        return "(use_spectrum)";
    }

    /// return the cosine of the angle between the incoming direction and the earth's zenith
    virtual double zenithCosine()const{return m_zenithCos;}

private:
    ISpectrum* m_spectrum;
    bool   m_galactic,m_equatorial;
    double m_zenithCos;
};
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//                   FluxSource constructor
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
FluxSource::FluxSource(const DOM_Element& xelem )
: EventSource ()
, m_spectrum(0)
, m_occultable(true)
, m_zenithCosTheta(1.0) //won't be occulted by default
{
    static double d2r = M_PI/180.;


    ISpectrum*   s = 0;
    std::string class_name;
    std::string source_params; 
    // this is a default flux, from the flux="123" in the source element
    setFlux(atof (xml::Dom::getAttribute(xelem, "flux").c_str()));


    DOM_Element   spec = xml::Dom::findFirstChildByName(xelem, "spectrum");

    if (spec == DOM_Element()) {

        // source has no imbedded spectrum element: expect a name
        class_name = xml::Dom::transToChar(xelem.getAttribute("name"));
    } else {
        // process spectrum element
        DOM_NodeList children = spec.getChildNodes();

        // First child element is type of spectrum
        DOM_Node    childNode = children.item(0);
        DOM_Element specType;

        if (childNode.getNodeType() == DOM_Node::ELEMENT_NODE) {
            specType = (DOM_Element &) childNode;
        }
        else specType = xml::Dom::getSiblingElement(childNode);

        DOMString   typeTagName = specType.getTagName();
        std::string spectrum_name = xml::Dom::transToChar(spec.getAttribute("name"));
        std::string spectrum_energyscale = xml::Dom::transToChar(spec.getAttribute("escale"));

        if(spectrum_energyscale == "GeV"){ m_energyscale=GeV;
        }else if(spectrum_energyscale == "MeV"){ m_energyscale=MeV;
        }else{
            std::cout << "bad energy scale declaration on spectrum:"
                << spectrum_energyscale << " , exiting.";
            return;} //this line "just in case"


        if (typeTagName.equals("particle")) s = new SimpleSpectrum(specType,  m_energyscale==GeV );
        else if (typeTagName.equals("SpectrumClass")) {
            // attribute "name" is the class name
            class_name = xml::Dom::transToChar(specType.getAttribute("name"));
            source_params= xml::Dom::transToChar(specType.getAttribute("params"));
        }
        else {
            // no, the tag itself
            class_name = xml::Dom::transToChar(typeTagName);//.transcode();
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
        }
        m_spectrum =s;


        // second child element is angle
        DOM_Element angles = xml::Dom::getSiblingElement(specType);
        DOMString anglesTag = angles.getTagName();

        if (anglesTag.equals("solid_angle") ) {
            m_occultable=false;
            m_launch_dir = new RandomDirection(
                atof(xml::Dom::transToChar(angles.getAttribute("mincos"))),
                atof(xml::Dom::transToChar(angles.getAttribute("maxcos"))),
                atof(xml::Dom::transToChar(angles.getAttribute("theta"))) * d2r, 
                atof(xml::Dom::transToChar(angles.getAttribute("phi"))) *d2r);

        }
        else if (anglesTag.equals("direction") ) {
            m_occultable=false;
            m_launch_dir = new LaunchDirection(
                atof(xml::Dom::transToChar(angles.getAttribute("theta"))) * d2r, 
                atof(xml::Dom::transToChar(angles.getAttribute("phi"))) *d2r);
        }
        else if (anglesTag.equals("use_spectrum") ) {
            std::string frame = xml::Dom::transToChar(angles.getAttribute("frame"));
            m_occultable=(frame=="galaxy" || frame=="equatorial");
            m_launch_dir = new SourceDirection(m_spectrum, frame); 
        }
        else if(anglesTag.equals("galactic_dir")){
            m_occultable=true;
            m_launch_dir = new LaunchDirection(
                astro::SkyDir(
                atof(xml::Dom::transToChar(angles.getAttribute("l"))),
                atof(xml::Dom::transToChar(angles.getAttribute("b"))), 
                astro::SkyDir::GALACTIC
                ),
                atof(xml::Dom::getAttribute(angles, "radius").c_str())
                );
        }
        else if(anglesTag.equals("celestial_dir")){
            m_occultable=true;
            m_launch_dir = new LaunchDirection(
                astro::SkyDir(
                atof(xml::Dom::transToChar(angles.getAttribute("ra"))),
                atof(xml::Dom::transToChar(angles.getAttribute("dec"))), 
                astro::SkyDir::CELESTIAL 
                ),

                atof(xml::Dom::getAttribute(angles, "radius").c_str())
                );

        }
        else if(anglesTag.equals("galactic_spread")){
            m_occultable=true;
            FATAL_MACRO("not implemented");
        }
        else {
            FATAL_MACRO("Unknown angle specification in Flux::Flux \""
                << xml::Dom::transToChar(anglesTag) << "\"" );
        }

        // third child element is optional launch spec
        DOM_Element launch = xml::Dom::getSiblingElement(angles);

        if(launch !=DOM_Element()) {
            DOMString launchTag = launch.getTagName();

            if(launchTag.equals("launch_point")){
                double float1=atof(xml::Dom::transToChar(launch.getAttribute("x")));
                double float2=atof(xml::Dom::transToChar(launch.getAttribute("y")));
                double float3=atof(xml::Dom::transToChar(launch.getAttribute("z")));
                m_launch_pt = new FixedPoint(HepPoint3D(float1,float2,float3),
                    atof(xml::Dom::transToChar(launch.getAttribute("beam_radius"))) );
            }else if(launchTag.equals("patch")){
                float num1=atof(xml::Dom::transToChar(launch.getAttribute("xmax")));
                float num2=atof(xml::Dom::transToChar(launch.getAttribute("xmin")));
                float num3=atof(xml::Dom::transToChar(launch.getAttribute("ymax")));
                float num4=atof(xml::Dom::transToChar(launch.getAttribute("ymin")));
                float num5=atof(xml::Dom::transToChar(launch.getAttribute("zmax")));
                float num6=atof(xml::Dom::transToChar(launch.getAttribute("zmin")));
                m_launch_pt = new Patch(num1,num2,num3,num4,num5,num6);
            }else {
                FATAL_MACRO("Unknown launch specification in Flux::Flux \""
                    << xml::Dom::transToChar(launchTag) << "\"" );
            }
        } else {
            // default: the target sphere.
            double radius = sqrt(totalArea() / M_PI ) * 1000;    // radius in mm
            m_launch_pt = new RandomPoint(radius, backoff_distance);
        }
    }
}


FluxSource::~FluxSource()
{
    delete m_spectrum;
    delete m_launch_pt;
    delete m_launch_dir;

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
    // Outputs - pointer to the "current" fluxSource object.
    m_interval = calculateInterval(time);
    computeLaunch(time+m_interval);
    //now set the actual interval to be what FluxMgr will get
    EventSource::setTime(time+m_interval);
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
    return -log(1.-RandFlat::shoot(1.))/rate(time);
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
    //HepRotation correctForTilt =GPS::instance()->rockingAngleTransform(time);
    m_correctedDir = /*correctForTilt*/m_launchDir;

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