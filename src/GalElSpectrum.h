// $Header: /nfs/slac/g/glast/ground/cvs/FluxSvc/src/GalElSpectrum.h,v 1.12 2003/02/23 02:08:22 burnett Exp $
// Original author: P. L. Nolan, pln@egret1.Stanford.EDU

/*
Notation: probably one of those things is called "flux" and the other
is something else.  Can someone look it up?
*/

#ifndef GAL_EL_SPECTRUM_H
#define GAL_EL_SPECTRUM_H

/** 
* \class GalElSpectrum
*
* \brief A quick and dirty implementation of the high-energy (galactic) cosmic ray spectrum.
* The spectrum is the power law found in Barwick et al., ApJ 498, 779 (1998).
* The low-energy rollover is ignored.
* Also ignored are the east-west effect and the gradual nature of the
* geomagnetic cutoff.  The power law is cut off sharply at the cutoff
* energy and directions are isotropic above the  horizon.
* The cutoff energy is obtained by querying a CHIMESpectrum object and
* adjusting for the smaller mass of the electron.  (Actually the electron
* mass is assumed to be zero.)

  *  As with CHIMESpectrum, the flux normalization is handled in a screwy
  *  way.  We should come up with a better way to handle this.
  *  The value returned is in electrons/(m^2 sec ster), as if it was
  *  isotropic.  This value is obtained by taking the observed isotropic
  *  flux and multiplying by a correction factor for the fraction of the
  *  sky blocked by the earth.  Thus the total integrated flux in
  *  electrons/(m^2 sec) can be obtained by multiplying by 4*pi.
  * \author P. L. Nolan, April 1999
  * 
  * $Header $
*/


#include "flux/Spectrum.h"
#include "facilities/Observer.h"
#include <string>
#include "CHIMESpectrum.h"

class GalElSpectrum : public Spectrum
{
    
public:
    /// ctor
    GalElSpectrum(const std::string& params);
    virtual ~GalElSpectrum();
    
    virtual double calculate_rate(double old_rate);
    
    /// calculate flux for the current cutoff
    virtual double flux(double) const;
    
    /// calcualte effective solid angle for the given energy 
    virtual double solidAngle()const;
    
    /// Total flux in electrons / m^2 sec ster
    float flux(float cut) const;
    
    virtual float flux(float lat, float lon) const;
    virtual float flux(std::pair<double, double> coords) const;
    
    virtual float operator() (float)const;
    
    float cutoff () const {return m_cutoff;};
    
    virtual void setPosition(float lat, float lon);
    virtual void setPosition(std::pair<double,double> coords);
    
    /// this one asks the GPS for position
    int askGPS();
    
    /// determine the cutoff value which will produce the desired flux
    float findCutoff(float rflux) const;
    
    /// determine the cutoff value at a geographical location
    float findCutoff(float lat, float lon) const;
    float findCutoff(std::pair<double,double> coords) const;
    
    /// return solid angle pair (costh, phi) for the given energy
    virtual std::pair<double,double> dir(double energy);
    
    
    virtual std::string title() const;
    virtual const char * particleName() const;
    inline  const char * nameOf() const {return "GalElSpectrum";}
    //   use default destructor, copy constructor, and assignment op.
    
    /// set the particle name. (default is "p")
    void setParticleName(std::string name);
    
private:
    CHIMESpectrum m_pspec;
    void init(float lat, float lon);
    float m_expo;   // exponent of power-law spectrum (integral)
    float m_norm;
    float m_normfact; // fraction of sky not blocked by earth
    float m_cutoff; // current cutoff energy
    float m_coscutoff;  // zenith angle of horizon
    float m_flux;   // current flux (set when cutoff changes)
    float m_lat, m_lon;  //current lat, lon.
    
    static const float m_rearth;    // radius of earth in km
    static const float m_altitude;  // altitude of circular orbit

    bool m_allowMove;   //whether or not the flux will be allowed to change position.
    
    std::string m_particle_name;
    ObserverAdapter< GalElSpectrum > m_observer; //obsever tag
};

#endif // GAL_EL_SPECTRUM_H


