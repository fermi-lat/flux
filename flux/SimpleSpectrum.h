/** @file SimpleSpectrum.h
    @brief declaration of SimpleSpectrum

   $Header: /nfs/slac/g/glast/ground/cvs/FluxSvc/src/SimpleSpectrum.h,v 1.6 2003/04/03 19:44:50 burnett Exp $
*/
#ifndef SIMPLESPECTRUM_H
#define SIMPLESPECTRUM_H

#include "Spectrum.h"
#include <string>

class DOM_Element;


/** 
* \class SimpleSpectrum
* @brief define a particle and spectral index
* 
* $Header: /nfs/slac/g/glast/ground/cvs/FluxSvc/src/SimpleSpectrum.h,v 1.6 2003/04/03 19:44:50 burnett Exp $
*/
class SimpleSpectrum : public Spectrum {
public: 
    SimpleSpectrum(const char* name,float E0, float index=0.0);
    SimpleSpectrum(const char* name,float Emin, float Emax, float index);
    SimpleSpectrum(const DOM_Element& xelem, bool useGeV=true);
    SimpleSpectrum(const std::string& params);
    
    SimpleSpectrum();
    void setPosition ( float /*lat*/, float /*lon*/ ){}
    virtual double calculate_rate(double old_rate);
    virtual float  operator()(float f)const;
    virtual const char* particleName()const;
    virtual std::string title()const;
private:
    float parseParamList(std::string input, int index);
    float m_E0;		// energy base
    std::string m_name;	// particle name to generate ("P", "gamma", ...)
    float m_index;	// spectral index: <=1 is delta function at E0
    float m_emax;
    bool m_useGeV;  // true if using GeV units, MeV otherwise
};

#endif
