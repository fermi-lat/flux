/** @file SimpleSpectrum.cxx
    @brief definition of SimpleSpectrum

   $Header: /nfs/slac/g/glast/ground/cvs/flux/src/SimpleSpectrum.cxx,v 1.4 2003/10/29 14:05:49 burnett Exp $
*/


#include "flux/SimpleSpectrum.h"

#include <xercesc/dom/DOM_Element.hpp>
#include "xml/Dom.h"

#include "flux/FluxException.h" // for FATAL_MACRO
#include <utility>
#include <sstream>
#include <cmath>
#include "flux/SpectrumFactory.h"

static SpectrumFactory<SimpleSpectrum> factory;
namespace {
    // useful utility functions

    // differential rate: return energy distrbuted as e**-gamma between e1 and e2, if r is uniform from 0 to1
    double power_law( double r, double e1, double e2, double gamma)
    {
       double e= gamma==1
           ?  e1*exp(r*log(e2/e1))
           :  e1*exp(log(1.0 - r*(1.-pow(e2/e1,1-gamma)))/(1-gamma));
        return e;
    }
    // integral of e**-gamma from e1 to e2
    double total_rate(double e1, double e2, double gamma)
    {
        return gamma==1
            ?  log(e2/e1)
            : ( pow(e1, 1-gamma) - pow(e2,1-gamma) ) / (gamma-1);
    }
}
    
// Is this used? if so, should be extended to set all parameters
SimpleSpectrum::SimpleSpectrum(const std::string& params)
:m_name("gamma")
,m_E0(parseParamList(params,0))
,m_index(parseParamList(params,1))
{}


SimpleSpectrum::SimpleSpectrum(const DOM_Element& xelem, bool useGeV)
: m_useGeV(useGeV)
{
    m_name = xml::Dom::getAttribute(xelem, "name").c_str();
    
    const DOM_Element spectrum = xml::Dom::findFirstChildByName(xelem, "*");
    
#if 0
    if (spectrum.getTagName().equals(DOMString("power_law"))) {
        m_E0 = atof(xml::Dom::getAttribute(spectrum, "emin").c_str());
        m_emax = atof(xml::Dom::getAttribute(spectrum, "emax").c_str());
        m_index = atof(xml::Dom::getAttribute(spectrum, "gamma").c_str());
        m_ebreak = atof(xml::Dom::getAttribute(spectrum, "ebreak").c_str());
        m_index2 =atof(xml::Dom::getAttribute(spectrum, "gamma2").c_str());
#else
    std::string tagName = xml::Dom::getTagName(spectrum);
    if( tagName == "power_law" ){
        m_E0 =      xml::Dom::getDoubleAttribute(spectrum, "emin");
        m_emax = xml::Dom::getDoubleAttribute(spectrum, "emax");
        m_index =xml::Dom::getDoubleAttribute(spectrum, "gamma"); 
        m_ebreak = xml::Dom::getDoubleAttribute(spectrum, "ebreak");
        m_index2 =xml::Dom::getDoubleAttribute(spectrum, "gamma2");

#endif
        if( m_ebreak==0) {
            m_ebreak=m_emax;
            m_a = 1.0; 
        }else{
            // calculate relative part of spectrum for lower
            double a1 = total_rate(m_E0, m_ebreak, m_index);
            double a2 =  pow(m_ebreak, m_index2-m_index)*total_rate(m_ebreak, m_emax, m_index2);
            m_a = a1/(a1+a2);
        }
    }
#if 0
    else if (spectrum.getTagName().equals(DOMString("energy"))) {
        m_E0 = atof(xml::Dom::getAttribute(spectrum, "e").c_str());
#else
    else if(tagName=="energy") {
#endif
        m_emax = 100.0;
        m_index = 0.0;
    }
#if 0
    else if (spectrum.getTagName().equals(DOMString("exponential"))) {
        m_E0 = atof(xml::Dom::getAttribute(spectrum, "emin").c_str());
        m_index = atof(xml::Dom::getAttribute(spectrum, "exponent").c_str());
#else
    else if( tagName == "exponential") {
        m_E0 = xml::Dom::getDoubleAttribute(spectrum, "exponential");
        m_index = xml::Dom::getDoubleAttribute(spectrum,"exponent");
#endif
        m_emax = 100.0;
        m_index = 0.0;
        FATAL_MACRO("exponential spectral component not implemented yet");
    }
    else {
        std::cerr << "Unknown name: " << m_name << std::endl;
        FATAL_MACRO("Unknown particle spectrum!");
    }
}


std::string SimpleSpectrum::title()const
{
    std::stringstream s;
    s << particleName() << '(' << m_E0 <<  (m_useGeV? " GeV" : " MeV") ;
    if( m_index >=1 ) s << ',' << m_index ;
    if(m_ebreak !=0) s << "," << m_ebreak << "," << m_index2;
    s << ")";
    return s.str();
}


float
SimpleSpectrum::operator()(float f)const
{
    if( m_index == 0.0 )     return m_E0;
    
    float energy;
    if( f<m_a ) {
        //float x = 1 - exp((1-m_index)*log(m_emax/m_E0));
        //return m_E0*exp(log(1-x*f)/(1-m_index));
        // single power law, or first segment
        energy =  power_law(f/m_a, m_E0, m_ebreak, m_index);
    }else{
        // break in the power law above the break
        energy = power_law( (f-m_a)/(1-m_a), m_ebreak, m_emax, m_index2);
    }
    return energy;
}

const char*
SimpleSpectrum::particleName() const
{
    return m_name.c_str();
}

float SimpleSpectrum::parseParamList(std::string input, int index)
{
    std::vector<float> output;
    int i=0;
    for(;!input.empty() && i!=std::string::npos;){
        float f = ::atof( input.c_str() );
        output.push_back(f);
        i=input.find_first_of(",");
        input= input.substr(i+1);
    } 
    return output[index];
}
