/** @file SimpleSpectrum.h
    @brief declaration of SimpleSpectrum

   $Header: /nfs/slac/g/glast/ground/cvs/flux/flux/SimpleSpectrum.h,v 1.3 2003/10/29 13:30:12 burnett Exp $
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
* $Header: /nfs/slac/g/glast/ground/cvs/flux/flux/SimpleSpectrum.h,v 1.3 2003/10/29 13:30:12 burnett Exp $
*/
class SimpleSpectrum : public Spectrum {
public: 

    /// ctor for instantiation from XML element  "particle"
    /// @param xelem nested element, expect either "power_law or "energy"
    SimpleSpectrum(const DOM_Element& xelem, bool useGeV=true);

    /// ctor for instantiation from XML element "SpectrumClass"
    /// @param parameter string 
    /// This is required, but the implementation has never been used
    SimpleSpectrum(const std::string& params);
    
    virtual float  operator()(float f)const;
    virtual const char* particleName()const;
    virtual std::string title()const;

    // convenient access methods
    double ebase()const{return m_E0;}
    double index()const{return m_index;}
    double ebreak()const{return m_ebreak;}
    double index2()const{return m_index2;}

private:
    float parseParamList(std::string input, int index);
    float m_E0;		// energy base
    std::string m_name;	// particle name to generate ("P", "gamma", ...)
    float m_index;	// spectral index: <=1 is delta function at E0
    float m_emax;
    float m_ebreak;    // if not zero, put in a break
    float m_index2;   // index for break
    bool m_useGeV;  // true if using GeV units, MeV otherwise
    double m_a;  // relative area of lower part of broken spectrum
};

#endif
