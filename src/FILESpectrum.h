/** 
* @file FILESpectrum.h
* @brief definition of FILESpectrum
*
*  $Header: /nfs/slac/g/glast/ground/cvs/FluxSvc/src/FILESpectrum.h,v 1.8 2002/11/20 20:09:54 srobinsn Exp $
*/
#ifndef FILESpectrum_H
#define FILESpectrum_H
/** 
* \class FILESpectrum
*
* \brief Spectrum that reads its differential spectrum from a table
* \author Theodore Hierath
* 
* $Header $:
*/
//
#include "flux/Spectrum.h"
#include "facilities/Observer.h"
#include <vector>
#include <string>
#include <map>


class FILESpectrum : public Spectrum
{
public:
    /// params is the filename to read
    FILESpectrum(const std::string& params);
    
    /// return total flux 
    virtual double flux() const;

    /// flux method to conform to FluxSvc spectrum standard
    double flux (double time ) const;
    
    /// sample a single particle energy from the spectrum
    virtual float operator() (float)const;
    
    virtual std::string title() const;
    virtual const char * particleName() const;
    inline  const char * nameOf() const {return "FILESpectrum";}
    //   use default destructor, copy constructor, and assignment op.
    
    
private:
    
    
    float m_flux;   // current flux (set when cutoff changes)
    
    
    typedef std::pair<double,double> efpair;
    std::vector<efpair> integ_flux;
    
    
    std::string initialization_document;
    std::string m_particle_name;
};



#endif // FILESpectrum_H
