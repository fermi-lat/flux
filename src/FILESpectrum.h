/** 
* @file FILESpectrum.h
* @brief FILESpectrum allows definition of a source based on a spectrum 
* saved in a file. The file is currently assumed to have two columns: energy value
* and flux (normalized or not). Flux value can be renormalized using the xml source 
* description. Energy in the file can be in MeV or GeV, also to be decided in the xml
* description.
*
*  $Header: /nfs/slac/g/glast/ground/cvs/users/cohen/flux/src/FILESpectrum.h,v 1.2 2004/09/27 21:43:10 cohen Exp $
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
    
    std::vector<std::pair<double,double> > initialize(const std::string& params);
    std::vector<std::pair<double,double> > readFile(std::ifstream& input_file);

private:
    
    
    float m_fileflux;   // normalization from the file: can be overwritten in the xml
    
    bool m_inLog;

    typedef std::pair<double,double> efpair;
    std::vector<efpair> integ_flux;
    
    
    std::string initialization_document;
};



#endif // FILESpectrum_H
