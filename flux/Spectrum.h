/** @file Spectrum.h
    @brief declaration of Spectrum
 
   $Header: /nfs/slac/g/glast/ground/cvs/FluxSvc/src/Spectrum.h,v 1.12 2003/03/02 06:30:57 burnett Exp $
*/

#ifndef GLAST_SPECTRUM_H
#define GLAST_SPECTRUM_H


#include <string>
#include <utility>
#include <vector>
#include "ISpectrum.h"


/** 
* \class Spectrum
* \brief base class for energy spectrum objects
*
* Spectrum is the base class for all of the particle sources 
* internal to FluxSvc.
 Class for holding function definitions of Spectrums - i.e. HeSpectrum, SimpleSpectrum, etc...
 Basically an abstract base class for these classes.
* 
* $Header: /nfs/slac/g/glast/ground/cvs/FluxSvc/src/Spectrum.h,v 1.12 2003/03/02 06:30:57 burnett Exp $
*/
class Spectrum : public ISpectrum {
public:
    
    //    class  Direction : public std::pair<float,float> {
    //    public:
    //        double costh()const{return this.first;}
    //        double phi()const {return this.second;}
    //    };
    
    virtual float operator()(float /*r*/)const{return 0;};
    // this is all that these objects do. Must be virtual for
    // polymorphism
    // returns kinetic energy for random number r in [0,1). 
    // requried that it be monatonic
    // NB. Default is to return zero, an indicator that the actual Spectrum object
    //     implements a method that makes direct use of the random generator
    
    //virtual double calculate_rate (double old_rate) = 0;
    
    /// subclasses need to specify correct particle type
    virtual const char * particleName()const=0;
    
    /// calculate the flux, particles/m^2/sr. (default zero)
    virtual double    flux (double time ) const;
    
    
    /// calcualte effective solid angle  (default zero)
    virtual double solidAngle()const;
    
    /// return a title describing the spectrum	
    virtual std::string title()const=0;
    
    /// inverse of the operator: given a KE, return number in [0,1)
    /// for choosing limits
    float fraction(float energy);
    
    virtual ~Spectrum();
    
    /// a randomized interval to the next event - default is 1/rate()
    virtual double interval (double time);
    
  
    /// new interface for Hirosima classes
    virtual double energy( double time=0);
    
    /*! 
    @param energy energy returned by previous call to energy
    \return dir direction is either in the format (cos theta, phi)
    (zenith-local coordinates, or (l,b) (galactic coordinates).
    */
    virtual std::pair<double,double> dir(double energy);
    
    
    
    
protected:
    Spectrum(const std::vector<float>& /*params*/){};
    Spectrum(){}
        // all constructors protected to ensure an abstract class
        
        virtual void parseParamList(std::string input, std::vector<float>& output) const;
    
    double m_currentInterval; // so we only find the interval for each particle once.
    
};

#endif    
