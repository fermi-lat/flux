/** @file FluxSource.h
@brief FluxSource declaration
*/
#ifndef FluxSource_h
#define FluxSource_h 1

#include "EventSource.h"
// forward declarations
#include "astro/SkyDir.h"

#include <xercesc/util/XercesDefs.hpp>
XERCES_CPP_NAMESPACE_BEGIN
class  DOMElement;
XERCES_CPP_NAMESPACE_END

class ISpectrum;
class LaunchDirection;
class LaunchPoint;

// 
/** @class FluxSource
@brief class which manages to compute flux from various particle source configurations
It is initialized from a xml description

$Header: /nfs/slac/g/glast/ground/cvs/flux/flux/FluxSource.h,v 1.11 2008/01/07 04:18:22 burnett Exp $
*/
class FluxSource : public EventSource  
{
public:      
    /**  constructor
    @param xelem The xml description for this source
    */
    FluxSource ( const XERCES_CPP_NAMESPACE_QUALIFIER DOMElement* xelem );

    ///    destructor
    virtual ~FluxSource();

    ///    generate an event 
    virtual EventSource* event(double time);
    ///    full-length title description of this EventSource.
    virtual std::string fullTitle () const;

    ///    brief title description (for display) for this event source
    virtual std::string displayTitle () const;

    virtual double flux(double time)const; // calculate flux for attached spectrum

    virtual double rate(double time)const; // calculate rate for attached spectrum

    /// return a title describing the spectrum and angles
    std::string title()const;

    /// print facility
    void  printOn ( std::ostream&  ) {}

    /// set spectrum, with optional parameter to set the maximum energy?
    void spectrum(ISpectrum* s, double emax=-1);

    ISpectrum* spectrum() const{ return m_spectrum; }

    

    //! Denotes what Energy Units the energy
    //! of incoming particles are in
    enum EnergyScale { 
        MeV,        //! MeV
        GeV         //! GeV
    } m_energyscale;

    virtual int eventNumber()const;

    virtual double energy()const { return m_energy;}
    virtual const CLHEP::Hep3Vector& launchDir()const {return m_correctedDir;}
    virtual const CLHEP::Hep3Vector&  launchPoint()const { return m_launchPoint;}
 
    virtual astro::SkyDir skyDirection()const;

    virtual std::string particleName();
    /// this function decides if the current incoming photon would be occulted
    /// by the earth
    bool occulted();

    virtual void disable();

private:

// forward declaration of nested classes that handle the lauch direction
   class RandomDirection;  // choose randomly from range 

// forward declaration of classes that handle launch point
   class RandomPoint; // random strategy
   class FixedPoint;  // fixed, or pencil
   class Patch;  // a box

    LaunchPoint* m_launch_pt; // pointer to actual point stategy: must be set
    LaunchDirection* m_launch_dir;

    ISpectrum*         m_spectrum;	    // spectrum to generate

    double m_energy;
    // associated with a specific launch

    /// result of strategy
    CLHEP::Hep3Vector m_launchDir;

    ///direction after being corrected for the "tilt" angles.
    CLHEP::Hep3Vector m_correctedDir;

    CLHEP::Hep3Vector  m_launchPoint;

    double calculateInterval (double time);

    ///interval function to be used by non-spectrum sources
    double explicitInterval (double time);
    ///    getLaunch - compute launch point, direction, & energy
    virtual void computeLaunch (double time=0);
    ///flag showing whether the current spectrum can be occulted by the earth.
    bool m_occultable;

    ///cosine of angle between zenith direction and incoming particle direction.
    double m_zenithCosTheta;

   bool m_launch_dir_owner;
   bool m_launch_pt_owner;

};
#endif
