// $Id: CHIMESpectrum.h,v 1.13 2003/02/23 02:08:22 burnett Exp $
#ifndef CHIME_SPECTRUM_H
#define CHIME_SPECTRUM_H

/** 
* \class CHIMESpectrum
*
* \brief Calculate the cosmic ray proton spectrum in low earth orbit.
*  Uses data produced by CHIME, assuming 600 km circular orbit
*  at solar minimum (worst case).
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
* Cutoff energies are returned by FindCutoff().  With no arguments,
* it returns the stored current value.  Two arguments are assumed
* to be latitude and longitude, and the corresponding cutoff is
* looked up.
* The operator() function returns a sampled energy value.  The argument
* must be a float value between 0 and 1.
* The dir() member returns a sampled particle's energy and direction.
* The argument is a value between 0 and 1, which is used to produce
* the energy, as in operator().  The direction is based on simple
* Stormer cone theory.  The first component of the direction is the
* zenith direction cosine: +1 for particles moving downward.  The 
* second component of the direction is the earth azimuth angle: 0
* from the east, +pi/2 from the north, -pi/2 from the south.  At
* energies near the cutoff, particles come mainly from the west.

  * \author Patrick Nolan, Stanford University, 1998
  * 
  * $Header $
*/

#include "flux/Spectrum.h"
#include "facilities/Observer.h"
#include <vector>
#include <string>


class CHIMESpectrum : public Spectrum
{
    
public:
    CHIMESpectrum(const std::string& params);
    
   
    /// calculate flux for the current cutoff
    virtual double flux(double) const;
    
    /// calcualte effective solid angle for the given energy 
    virtual double solidAngle()const;
    
    
    /// Total flux in protons / m^2 sec ster
    /// Interpolate in table if possible, otherwise assume power law
    ///  tail at high energy.
    float flux(float cut) const;
    
    /// Flux as a function of latitude and longitude in a 600 km orbit.
    /// Linear interpolate in a table with a 5 degree sampling grid.
    virtual float flux() const;
    
    virtual float operator() (float)const;
    
    float cutoff () const {return m_cutoff;};
    
    virtual void setPosition(double lat, double lon);
    virtual void setPosition(std::pair<double,double> coords);
    
    /// this one asks the GPS for position
    int askGPS();
    
    /// determine the cutoff value which will produce the desired flux
    float findCutoff(float rflux) const;
    
    /// determine the cutoff value at the current geographical location
    float findCutoff() const;
    
    /// return solid angle pair (costh, phi) for the given energy
    virtual std::pair<double,double> dir(double energy);
    
    virtual std::string title() const;
    virtual const char * particleName() const;
    inline  const char * nameOf() const {return "CHIMESpectrum";}
    //   use default destructor, copy constructor, and assignment op.
    
    /// set the particle name. (default is "p")
    void setParticleName(std::string name);
    
    
    typedef std::pair<int,float> Intrp;
    
protected:
    void init(std::string params);
private:
    /// a vector with methods added to do searches and linear interpolation
    class InterpVec : public std::vector<float> {
        
    public:
        InterpVec();  // constructor
        Intrp search(float x) const;
        float interpolate(Intrp) const;
    };
    
    
    InterpVec m_en; // values of energy in the table
    InterpVec m_fl; // values of integral flux in the table
    InterpVec m_fluxes; // integral flux, unmodulated by earth's field
    float m_expo;   // exponent of power-law spectrum above table (integral)
    double /*float*/ m_normfact; // fraction of sky not blocked by earth
    float m_tot;    // total flux above the current cutoff energy
    float m_upper;  // total flux in the power-law part of spectrum
    float m_etop;   // energy at which table transitions to power law
    float m_cutoff; // current cutoff energy
    float m_coscutoff;  // zenith angle of horizon
    float m_fluxTbl[73][13];  // table of total flux vs. latitude and longitude
    float m_flux;   // current flux (set when cutoff changes)
    
    float rad2() const; // square of distance (in km) from dipole center
    float exposure(float E) const;  // geomagnetic cutoff factor
    float cosomega(float E) const;  // Opening angle of Størmer cone
    static const float m_rearth;    // radius of earth in km
    static const float m_altitude;  // altitude of circular orbit
    bool m_allowMove;   //whether or not the flux will be allowed to change position.
    
    std::string m_particle_name;
    ObserverAdapter< CHIMESpectrum > m_observer; //obsever tag
    
    float m_lat, m_lon;
};

#endif // CHIME_SPECTRUM_H

