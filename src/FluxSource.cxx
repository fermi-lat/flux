/** @file FluxSource.cpp
@brief Implementation of FluxSource
*/

#include "astro/SkyDir.h"
#include "astro/GPS.h"

#include "flux/FluxSource.h"
#include "flux/LaunchDirection.h"
#include "flux/LaunchPoint.h"
#include "flux/SpectrumFactoryTable.h"
#include "flux/SimpleSpectrum.h"
#include "flux/FluxException.h" // for FATAL_MACRO

#include "SourceDirection.h"  // This defines SourceDirection class

// RapidXML includes
#include "xmlBase/rapidxml.hpp"

#include "CLHEP/Random/RandFlat.h"

#include <algorithm>
#include <sstream>
#include <stdexcept>
#include <cstring>
#include <cmath>
#include <list>

// Helper function implementations
std::string FluxSource::getAttribute(XmlNode* node, const char* name) {
    if (!node) return "";
    auto* attr = node->first_attribute(name);
    return attr ? std::string(attr->value()) : "";
}

double FluxSource::getDoubleAttribute(XmlNode* node, const char* name) {
    std::string val = getAttribute(node, name);
    if (val.empty()) return 0.0;
    try {
        return std::stod(val);
    } catch (...) {
        return 0.0;
    }
}

XmlNode* FluxSource::findFirstChildByName(XmlNode* parent, const char* name) {
    if (!parent) return nullptr;
    // Handle wildcard "*" to get first child element
    if (name && std::strcmp(name, "*") == 0) {
        return getFirstChildElement(parent);
    }
    return parent->first_node(name);
}

XmlNode* FluxSource::getFirstChildElement(XmlNode* parent) {
    if (!parent) return nullptr;
    for (auto* child = parent->first_node(); child; child = child->next_sibling()) {
        if (child->type() == rapidxml::node_element) {
            return child;
        }
    }
    return nullptr;
}

std::string FluxSource::getTagName(XmlNode* node) {
    if (!node || !node->name()) return "";
    return std::string(node->name(), node->name_size());
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
/** @class RandomPoint
@brief nested launch strategy derived class
This is the standard strategy, which takes a direction and creates a point in
a disk centered at the origin with area 6 m^2 (or so)
*/
class FluxSource::RandomPoint : public LaunchPoint { 
public:
    RandomPoint(double radius, double backoff)
        : m_radius(radius), m_backoff(backoff)
    { 
    }

    virtual void execute(const CLHEP::Hep3Vector& dir) {
        CLHEP::HepRotation r_pln;

        // create rotation to take x-y plane to be perpendicular to incoming direction
        double ly = dir.y(), lx = dir.x();
        if (fabs(lx) + fabs(ly) > 1e-8) {  // leave as identity 
            r_pln.rotate(acos(dir.z()), CLHEP::Hep3Vector(-ly, lx, 0.));
        }

        // pick a random position on the planar section of a sphere through 
        // its midpoint
        double azimuth = CLHEP::RandFlat::shoot(2 * M_PI);
        double rad = m_radius * sqrt(CLHEP::RandFlat::shoot());

        CLHEP::Hep3Vector posLaunch(rad * cos(azimuth), rad * sin(azimuth), 0.);

        // define actual launch point
        setPoint(r_pln * posLaunch - m_backoff * dir);
    }

    /// return info
    virtual std::string title() const {
        std::stringstream t;
        t << "radius(" << m_radius << ")";
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
class FluxSource::FixedPoint : public LaunchPoint { 
public:
    FixedPoint(const CLHEP::Hep3Vector& pt, double radius)
        : LaunchPoint(pt)
        , m_disk_radius(radius)
        , m_base_point(pt)
    {}

    virtual void execute(const CLHEP::Hep3Vector& dir) {
        if (m_disk_radius == 0) return; // just use base point

        CLHEP::HepRotation r_pln;

        double ly = dir.y(), lx = dir.x();
        if (lx != 0 || ly != 0) { 
            r_pln.rotate(acos(dir.z()), CLHEP::Hep3Vector(-ly, lx, 0.));
        }
        double azimuth = CLHEP::RandFlat::shoot(2 * M_PI);
        double rad = m_disk_radius * sqrt(CLHEP::RandFlat::shoot());
        CLHEP::Hep3Vector posLaunch(rad * cos(azimuth), rad * sin(azimuth), 0.);

        setPoint(r_pln * posLaunch + m_base_point);
    }

    virtual std::string title() const {
        if (m_disk_radius == 0) return LaunchPoint::title();
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
    Patch(double xmin, double xmax, double ymin, double ymax, double zmin, double zmax)
        : m_xmin(xmin), m_dx(xmax - xmin)
        , m_ymin(ymin), m_dy(ymax - ymin)
        , m_zmin(zmin), m_dz(zmax - zmin)
    {
    }

    virtual void execute(const CLHEP::Hep3Vector&) {
        setPoint(CLHEP::Hep3Vector(
            m_xmin + m_dx * CLHEP::RandFlat::shoot(),
            m_ymin + m_dy * CLHEP::RandFlat::shoot(),
            m_zmin + m_dz * CLHEP::RandFlat::shoot()));
    }

    virtual std::string title() const {
        std::stringstream t;
        t << "patch(" 
          << m_xmin << "," << m_xmin + m_dx << ","
          << m_ymin << "," << m_ymin + m_dy << ","
          << m_zmin << "," << m_zmin + m_dz << ")";
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
class FluxSource::RandomDirection : public LaunchDirection { 
public:
    /** ctor:
    @param minc  minimum value of cos(theta)
    @param maxc  maximum value of cos(theta)
    @param theta [0] X rotation angle (radians)
    @param phi   [0] Z rotation angle (radians)
    */
    RandomDirection(double minc, double maxc, double theta = 0, double phi = 0)
        : m_theta(theta)
        , m_phi(phi)
    {
        using std::min;
        using std::max;

        // require _maxCos > _minCos for solid angle calculation
        m_minCos = min(minc, maxc);
        m_maxCos = max(minc, maxc);
        if (m_minCos == m_maxCos) {
            if (m_minCos != -1) m_minCos -= 0.001; 
            else m_maxCos += 0.001;
        }
        m_minPhi = 0; 
        m_maxPhi = 2 * M_PI;
    }

    virtual void execute(double /*ke*/, double time) override {
        double costh = -CLHEP::RandFlat::shoot(m_minCos, m_maxCos);
        double sinth = sqrt(1. - costh * costh);
        double phi = CLHEP::RandFlat::shoot(m_minPhi, m_maxPhi);

        // here, the direction is with respect to the zenith frame,
        // so we need the transformation from the zenith to GLAST.
        CLHEP::HepRotation zenToGlast = astro::GPS::instance()->transformToGlast(time, astro::GPS::ZENITH);

        CLHEP::Hep3Vector dir(cos(phi) * sinth, sin(phi) * sinth, costh);

        // extra rotation in case not zenith pointing
        if (m_theta != 0.0) dir.rotateX(m_theta).rotateZ(m_phi);
        
        // Set the protected member m_lat_dir directly
        m_lat_dir = zenToGlast * dir;
    }

    //! solid angle
    virtual double solidAngle() const override {
        return 2 * M_PI * (m_maxCos - m_minCos);
    }

    virtual std::string title() const override {
        std::stringstream t;
        t << "range(" << m_minCos << ',' << m_maxCos << ") ";
        if (m_theta != 0) {
            t << ", angle(" << m_theta * 180 / M_PI << ',' << m_phi * 180 / M_PI << ") ";
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
FluxSource::FluxSource(XmlNode* xelem)
    : EventSource()
    , m_spectrum(0)
    , m_occultable(true)
    , m_zenithCosTheta(1.0)  // won't be occulted by default
    , m_launch_dir_owner(true)
    , m_launch_pt_owner(true)
{
    static double d2r = M_PI / 180.;
    setName(getAttribute(xelem, "name").c_str());

    ISpectrum* s = 0;
    std::string class_name;
    std::string source_params; 
    
    // this is a default flux, from the flux="123" in the source element
    setFlux(atof(getAttribute(xelem, "flux").c_str()));
    int ident(static_cast<int>(atof(getAttribute(xelem, "ident").c_str())));

    XmlNode* spec = findFirstChildByName(xelem, "spectrum");

    if (spec == 0) {
        // source has no embedded spectrum element: expect a name
        class_name = getAttribute(xelem, "name");
    } else {
        // process spectrum element
        XmlNode* specType = getFirstChildElement(spec);
        
        std::string typeTagName = getTagName(specType);
        //std::string particle_name = getAttribute(spec, "particle_name");
	XmlNode* particle = findFirstChildByName(spec, "particle");
	std::string particle_name = getAttribute(particle, "name");
	std::string spectrum_energyscale = getAttribute(spec, "escale");

        std::string apply_edisp = getAttribute(spec, "apply_edisp");
        if (apply_edisp != "true" && apply_edisp != "false" && apply_edisp != "") {
            throw std::runtime_error("Invalid value for apply_edisp attribute in xml definition of " + name());
        }
        m_applyEdisp = (apply_edisp != "false");

        if (spectrum_energyscale == "GeV") { 
            m_energyscale = GeV;
        } else if (spectrum_energyscale == "MeV") { 
            m_energyscale = MeV;
        } else {
            std::cout << "bad energy scale declaration on spectrum:"
                      << spectrum_energyscale << " , exiting.";
            return;
        }

        if (typeTagName == "particle") {
            s = new SimpleSpectrum(specType, m_energyscale == GeV);
        } else if (typeTagName == "SpectrumClass") {
            // attribute "name" is the class name
            class_name = getAttribute(specType, "name");
            source_params = getAttribute(specType, "params");
        } else {
            // no, the tag itself
            class_name = typeTagName;
        }

        // if s is still 0, we need to create the internal spectrum object.
        if (s == 0) {
            s = SpectrumFactoryTable::instance()->instantiate(class_name, source_params);
            if (s == 0) {
                std::cerr << "List of known Spectrum classes:\n";
                std::list<std::string> list = SpectrumFactoryTable::instance()->spectrumList();
                for (auto it = list.begin(); it != list.end(); ++it)
                    std::cerr << "\t" << *it << std::endl;
                FATAL_MACRO("Unknown Spectrum: " << class_name);
                return;
            }
            std::string flux = getAttribute(spec, "flux");
            if (!flux.empty()) {
                s->setFlux(atof(flux.c_str()));
            }
            s->setInGeV(spectrum_energyscale == "GeV");

            if (!particle_name.empty()) s->setParticleName(particle_name);
        }
        m_spectrum = s;
        m_spectrum->setIdentifier(ident);
    }

    // Process direction/angles element - these are children of the spectrum element
    XmlNode* angles = nullptr;
    if (spec) {
        // Look for angle specifications as children of spectrum (after SpectrumClass/particle)
        angles = findFirstChildByName(spec, "solid_angle");
        if (!angles) angles = findFirstChildByName(spec, "direction");
        if (!angles) angles = findFirstChildByName(spec, "use_spectrum");
        if (!angles) angles = findFirstChildByName(spec, "galactic_dir");
        if (!angles) angles = findFirstChildByName(spec, "celestial_dir");
        if (!angles) angles = findFirstChildByName(spec, "custom_dir");
    }

    if (angles) {
        std::string anglesTag = getTagName(angles);
        
        if (anglesTag == "solid_angle") {
            m_occultable = false;
            m_launch_dir = new RandomDirection(
                getDoubleAttribute(angles, "mincos"),
                getDoubleAttribute(angles, "maxcos"),
                getDoubleAttribute(angles, "theta") * d2r,
                getDoubleAttribute(angles, "phi") * d2r);
        }
        else if (anglesTag == "direction") {
            std::string frame = getAttribute(angles, "frame");
            m_occultable = (frame == "zenith");
            m_launch_dir = new LaunchDirection(
                getDoubleAttribute(angles, "theta") * d2r,
                getDoubleAttribute(angles, "phi") * d2r,
                frame);
        }
        else if (anglesTag == "use_spectrum") {
            std::string frame = getAttribute(angles, "frame");
            m_occultable = (frame != "zenith");
	    if (frame.empty()) frame = "zenith";  // DEFAULT IS SET HERE
            m_launch_dir = new SourceDirection(m_spectrum, frame);
        }
        else if (anglesTag == "galactic_dir") {
            m_occultable = true;
            m_launch_dir = new LaunchDirection(
                astro::SkyDir(
                    getDoubleAttribute(angles, "l"),
                    getDoubleAttribute(angles, "b"),
                    astro::SkyDir::GALACTIC),
                getDoubleAttribute(angles, "radius"));
        }
        else if (anglesTag == "celestial_dir") {
            m_occultable = true;
            m_launch_dir = new LaunchDirection(
                astro::SkyDir(
                    getDoubleAttribute(angles, "ra"),
                    getDoubleAttribute(angles, "dec")),
                getDoubleAttribute(angles, "radius"));
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
                        << anglesTag << "\"");
        }
    } else {
        // Default: random from hemisphere
        m_launch_dir = new RandomDirection(-1.0, 1.0);
        m_occultable = false;
    }

    // Process launch point element - look for third child after angles
    XmlNode* launch = nullptr;
    if (angles) {
        for (auto* sibling = angles->next_sibling(); sibling; sibling = sibling->next_sibling()) {
            if (sibling->type() == rapidxml::node_element) {
                launch = sibling;
                break;
            }
        }
    }
    if (!launch) {
        launch = findFirstChildByName(xelem, "launch_point");
        if (!launch) launch = findFirstChildByName(xelem, "patch");
        if (!launch) launch = findFirstChildByName(xelem, "custom_pt");
    }

    if (launch) {
        std::string launchTag = getTagName(launch);
        
        if (launchTag == "launch_point") {
            m_launch_pt = new FixedPoint(CLHEP::Hep3Vector(
                getDoubleAttribute(launch, "x"),
                getDoubleAttribute(launch, "y"),
                getDoubleAttribute(launch, "z")),
                getDoubleAttribute(launch, "beam_radius"));
        }
        else if (launchTag == "patch") {
            m_launch_pt = new Patch(
                getDoubleAttribute(launch, "xmax"),
                getDoubleAttribute(launch, "xmin"),
                getDoubleAttribute(launch, "ymax"),
                getDoubleAttribute(launch, "ymin"),
                getDoubleAttribute(launch, "zmax"),
                getDoubleAttribute(launch, "zmin"));
        }
        else if (launchTag == "custom_pt") {
            m_launch_pt = m_spectrum->launchPoint();
            if (m_launch_pt == 0) {
                std::ostringstream what;
                what << "FluxSource: cannot use a 'custom_pt' tag with a "
                     << class_name << " source.";
                throw std::runtime_error(what.str());
            }
            m_launch_pt_owner = false;
        }
        else {
            FATAL_MACRO("Unknown launch specification in Flux::Flux \""
                        << launchTag << "\"");
        }
    } else {
        // Default: random point on target sphere
        double radius = sqrt(totalArea() / M_PI) * 1000;   // radius in mm
        m_launch_pt = new RandomPoint(radius, EventSource::s_backoff);
    }
}

FluxSource::~FluxSource() {
    delete m_spectrum;
    if (m_launch_pt_owner) delete m_launch_pt;
    if (m_launch_dir_owner) delete m_launch_dir;
}

void FluxSource::spectrum(ISpectrum* s, double emax) {
    if (emax > 0) {
        std::cerr << "exercising obsolete function fraction" << std::endl;
    }
    m_spectrum = s;
}

EventSource* FluxSource::event(double time) {
    // Purpose and Method: generate a new incoming particle
    // Inputs  - current time
    // Outputs - pointer to the "current" fluxSource object, or zero if it has "turned off"
    if (!enabled()) {
        throw std::runtime_error("FluxSource::event called when disabled");
    }
    using astro::GPS;
    setInterval(calculateInterval(time));
    if (interval() <= 0) {
        throw std::runtime_error("EventSource::event: negative or zero interval");
    }
    if (time + interval() < GPS::instance()->endTime()) {
        // do this only if in valid interval: assume will never get used otherwise
        computeLaunch(time + interval());
    } else {
        // flag to end use of this source
        disable();
    }
    return this;
}

double FluxSource::calculateInterval(double time) {
    if (m_spectrum) {
        return m_spectrum->interval(time);
    }
    return explicitInterval(time);
}

double FluxSource::explicitInterval(double time) {
    // Default implementation - exponential interval based on rate
    double r = rate(time);
    if (r <= 0) return 1e10;  // Very large interval
    return -log(CLHEP::RandFlat::shoot()) / r;
}

void FluxSource::computeLaunch(double time) {
    if (m_spectrum) {
        // Use the energy() method from ISpectrum interface
        m_energy = m_spectrum->energy(time);
    }
    
    m_launch_dir->execute(m_energy, time);
    m_correctedDir = m_launch_dir->dir();
    
    // Get zenith cos theta from launch direction using zenithCosine() method
    m_zenithCosTheta = m_launch_dir->zenithCosine();
    
    m_launch_pt->execute(m_correctedDir);
    m_launchPoint = (*m_launch_pt)();
}

std::string FluxSource::fullTitle() const {
    return title();
}

std::string FluxSource::displayTitle() const {
    std::stringstream s;
    s << EventSource::displayTitle() << '(' << m_spectrum->title();
    s << ')' << '\0';
    return s.str();
}

int FluxSource::eventNumber() const {
    return 0;
}

std::string FluxSource::title() const {
    if (m_spectrum == 0) return "";
    return m_spectrum->title() + ", "
        + m_launch_pt->title() + ", "
        + m_launch_dir->title();
}

std::string FluxSource::particleName() {
    return spectrum()->particleName();
}

bool FluxSource::occulted() {
    using astro::GPS;
    using astro::SkyDir;
    // Purpose: to determine whether or not the current incoming particle will be blocked by the earth.
    
    static double minCosTheta = -0.4;

    if (!m_occultable) return false;
    if (m_zenithCosTheta < minCosTheta) return true;
    if (EventSource::s_cone.size() < 3) return false;
    
    // Additional cone filtering logic would go here if needed
    return false;
}

double FluxSource::flux(double time) const {
    if (m_spectrum) {
        return m_spectrum->flux(time);
    }
    return EventSource::flux(time);
}

double FluxSource::rate(double time) const {
    // Calculate rate from flux, solid angle, and total area
    double solidAng = m_launch_dir ? m_launch_dir->solidAngle() : 1.0;
    return flux(time) * solidAng * totalArea();
}

astro::SkyDir FluxSource::skyDirection() const {
    // Convert the LAT direction back to a SkyDir
    // The dir() method returns direction in LAT frame
    return astro::GPS::instance()->toSky(m_launch_dir->dir());
}

void FluxSource::disable() {
    m_enabled = false;
}

int FluxSource::identifier() {
    if (m_spectrum) {
        return m_spectrum->identifier();
    }
    return -1;
}

std::string FluxSource::name() const {
    return EventSource::name();
}
