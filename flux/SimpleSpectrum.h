/** @file SimpleSpectrum.h
    @brief declaration of SimpleSpectrum
*/
#ifndef SIMPLESPECTRUM_H
#define SIMPLESPECTRUM_H

#include "Spectrum.h"
#include <string>

// RapidXML forward declaration
namespace rapidxml {
    template<class Ch> class xml_node;
}

// Type alias for RapidXML node
using XmlNode = rapidxml::xml_node<char>;

/** 
* \class SimpleSpectrum
* @brief define a particle and spectral index
*/
class SimpleSpectrum : public Spectrum {
public: 

    /// ctor for instantiation from XML element "particle"
    /// @param xelem nested element, expect either "power_law" or "energy"
    SimpleSpectrum(XmlNode* xelem, bool useGeV = true);

    /// ctor for instantiation from XML element "SpectrumClass"
    /// @param parameter string 
    SimpleSpectrum(const std::string& params);
    
    virtual float operator()(float f);
    virtual const char* particleName() const;
    virtual std::string title() const;

    // convenient access methods
    double ebase() const { return m_E0; }
    double index() const { return m_index; }
    double ebreak() const { return m_ebreak; }
    double index2() const { return m_index2; }

private:
    float parseParamList(std::string input, int index);
    void setup_power_law();
    
    float m_E0;         // energy base
    std::string m_name; // particle name to generate ("P", "gamma", ...)
    float m_index;      // spectral index: <=1 is delta function at E0
    float m_emax;
    float m_ebreak;     // if not zero, put in a break
    float m_index2;     // index for break
    bool m_useGeV;      // true if using GeV units, MeV otherwise
    double m_a;         // relative area of lower part of broken spectrum
    
    // Helper functions for XML operations
    static std::string getAttribute(XmlNode* node, const char* name);
    static double getDoubleAttribute(XmlNode* node, const char* name);
    static XmlNode* findFirstChildByName(XmlNode* parent, const char* name);
    static std::string getTagName(XmlNode* node);
};

#endif
