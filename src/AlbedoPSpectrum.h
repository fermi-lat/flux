// $Id: AlbedoPSpectrum.h,v 1.11 2003/02/23 02:08:22 burnett Exp $
// File: AlbedoPSpectrum.h
#ifndef ALBEDO_P_SPECTRUM_H
#define ALBEDO_P_SPECTRUM_H

/** 
* \class AlbedoPSpectrum
*
* \brief Calculate the earth albedo proton spectrum in low earth orbit.
* Uses data produced by AMS detector, preliminary graphs from web
* page.  No angular or geographic dependence included.  
*
* Interface:
* The constructor arguments specify the satellite position
* (latitude, longitude) in degrees.
* The position can be changed by use of the setPosition() member.
* The total flux in protons/(m^2 sec ster) is returned by the
* flux() member.  There are 3 ways to call flux(): With no arguments
* it uses the current cutoff energy.  If there is one argument, it
* is assumed to be a cutoff.  Two arguments are assumed to be latitude
* and longitude, and the corresponding cutoff is looked up.  The
* stored current value is not changed.
* The operator() function returns a sampled energy value.  The argument
* must be a float value between 0 and 1.
* The dir() member returns a sampled particle's direction.  The 
* argument is the particle energy, as produced by ().

  * \author Patrick Nolan, Stanford University, 1999
  * 
  * $Header $
*/

#include "flux/Spectrum.h"
#include "facilities/Observer.h"
#include <string>


class AlbedoPSpectrum : public Spectrum
{
    
public:
    AlbedoPSpectrum(const std::string& params);
    
    virtual double calculate_rate(double old_rate);
    
    /// calculate flux for the current position
    virtual double flux(double) const;
    
    /// effective solid angle for the given energy
    virtual double solidAngle()const;
    
    /// flux as a function of latitude and longitude in 600 km orbit
    virtual float flux(float lat, float lon) const;
    virtual float flux(std::pair<double, double> coords) const;
    
    /// sample a single particle energy from the spectrum
    virtual float operator() (float)const;
    
    /// move to a new position and do the necessary initialization
    virtual void setPosition(float lat, float lon);
    virtual void setPosition(std::pair<double,double> coords);
    int askGPS(); // this one asks the GPS for position
    
    /// sample a solid angle pair (costh,phi) from angular distribution
    virtual std::pair<double,double> dir(double energy);
    
    virtual std::string title() const;
    virtual const char * particleName() const;
    inline  const char * nameOf() const {return "AlbedoPSpectrum";}
    //   use default destructor, copy constructor, and assignment op.
    
    /// set the particle name. (default is "p")
    void setParticleName(std::string name);
    
    
private:
    
    void init(const std::vector<float>& params);
    
    /// current flux (set when cutoff changes)
    float m_flux;   

    double m_lat,m_lon; // latitude, longitude
    bool m_allowMove;   //whether or not the flux will be allowed to change position.

    std::string m_particle_name;
    ObserverAdapter< AlbedoPSpectrum > m_observer; //obsever tag
    void fitParams(const double lat, const double lon, double& alf1, 
        double& alf2, double& emin, double& emax, double& v1, double& v2, 
        double& Ejoin) const;
};

#endif // ALBEDO_P_SPECTRUM_H

