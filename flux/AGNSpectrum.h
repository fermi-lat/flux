// $Header: /nfs/slac/g/glast/ground/cvs/flux/flux/AGNSpectrum.h,v 1.2 2004/09/24 22:44:15 burnett Exp $


#ifndef AGN_SPECTRUM_H
#define AGN_SPECTRUM_H

/** 
* \class AGNSpectrum
*
* \brief base class for AGN spectrum objects
*
* AGNSpectrum represents AGN sources, given input nominal flux, index,  
* flaring period, flaring flux multiplier and spectral index.
* 
* $Header: /nfs/slac/g/glast/ground/cvs/flux/flux/AGNSpectrum.h,v 1.2 2004/09/24 22:44:15 burnett Exp $
*/

#include <string>
#include <utility>
#include <vector>
#include "Spectrum.h"
#include <cmath>

#include "CLHEP/Random/RandFlat.h"
//forward declaration
class HepRandomEngine;


//!  Class for holding function definitions of Spectrums - i.e. HeSpectrum, SimpleSpectrum, etc...
//!  Basically an abstract base class for these classes.
class AGNSpectrum : public Spectrum {
public:
    
    /// Return energy
    virtual float operator()(float r);
    
    /// Constructor with parameters
    AGNSpectrum (const std::string& params);
    
    /// constructor with individual parameters
    AGNSpectrum (float l = 0.0f,float b = 0.0f,float flux = 10.0f,float index = 2.15f,float flareMult = 5.0f,
        float flareAdd = 0.01f,float flareDuty = 0.00f,float flarePeriod = 1.0f,float Emin=0.1f,float Emax=100.0f);
    
    /// subclasses need to specify correct particle type
    virtual const char * particleName()const{return "gamma";}
    
    /// calculate the flux, particles/m^2/sr. (default zero)
    virtual double    flux (double time ) const;
    
    /// calcualte effective solid angle  (default zero)
    virtual double solidAngle()const;
    
    /// return a title describing the spectrum	
    virtual std::string title()const{return "AGN";}
    
    virtual ~AGNSpectrum();
    
    /// a randomized interval to the next event - default is 1/rate()
    virtual double interval (double time);
    
    /// new interface for Hirosima classes
    virtual double energy(double time=0);
    
    /*! 
    \param dir direction is either in the format (cos theta, phi)
    (zenith-local coordinates, or (l,b) (galactic coordinates).
    */
    virtual std::pair<double,double> dir(double energy);
    
    ///function for handling initialization through the xml "params" variable.
    void init(std::vector<float> params);
    
    
    ///turns flaring to the opposite state.
    void flipState(){
        if(m_flaring){m_flaring=false;
        }else{m_flaring=true;
        }
    }
    
    /// set the initial flaring state of the AGN, and when it will flare next.
    void setInitialFlareStates();
    
    /// recalculate the next flare state flip time.
    void reCalcNextFlip();
    
protected:
    
    virtual void parseParamList(std::string input, std::vector<float>& output) const;
    
    float m_l, m_b;
    float m_flux;
    float m_index;
    float m_flareMult;
    float m_flareAdd;
    float m_flareDuty;
    float m_flarePeriod;
    float m_Emin,m_Emax; //energy extrema
    
    bool m_flaring; //is this AGN flaring?        
    double m_nextFlipTime;//the next time when the source flares or goes quiescent.
    
};

#endif    
