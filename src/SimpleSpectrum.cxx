/** @file SimpleSpectrum.cxx
    @brief definition of SimpleSpectrum

   $Header: /nfs/slac/g/glast/ground/cvs/FluxSvc/src/SimpleSpectrum.cxx,v 1.11 2003/04/03 19:44:50 burnett Exp $
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


SimpleSpectrum::SimpleSpectrum(){}//default constructor
SimpleSpectrum::SimpleSpectrum(const std::string& params)
:m_name("gamma")
,m_E0(parseParamList(params,0))
,m_index(parseParamList(params,1))
{}


SimpleSpectrum::SimpleSpectrum(const char* name, float E0, float index)
:m_E0(E0)
,m_name(name)
,m_index( index)
,m_emax(100.)
,m_useGeV(true)
{}

SimpleSpectrum::SimpleSpectrum(const char* name, float Emin, float Emax, float index)
:m_E0(Emin)
,m_name(name)
,m_index(index)
,m_emax(Emax)
,m_useGeV(true)
{}

SimpleSpectrum::SimpleSpectrum(const DOM_Element& xelem, bool useGeV)
: m_useGeV(useGeV)
{
    m_name = xml::Dom::getAttribute(xelem, "name").c_str();
    
    const DOM_Element spectrum = xml::Dom::findFirstChildByName(xelem, "*");
    if (spectrum.getTagName().equals(DOMString("power_law"))) {
        m_E0 = atof(xml::Dom::getAttribute(spectrum, "emin").c_str());
        m_emax = atof(xml::Dom::getAttribute(spectrum, "emax").c_str());
        m_index = atof(xml::Dom::getAttribute(spectrum, "gamma").c_str());
    }
    else if (spectrum.getTagName().equals(DOMString("energy"))) {
        m_E0 = atof(xml::Dom::getAttribute(spectrum, "e").c_str());
        m_emax = 100.0;
        m_index = 0.0;
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
    s << ")";
    return s.str();
}

float
SimpleSpectrum::operator()(float f)const
{
    if( m_index == 0.0 )     return m_E0;
    
    if( m_index == 1.0 ) return m_E0*exp(f*log(m_emax/m_E0));
    
    float x = 1 - exp((1-m_index)*log(m_emax/m_E0));
    return m_E0*exp(log(1-x*f)/(1-m_index));
}

const char*
SimpleSpectrum::particleName() const
{
    return m_name.c_str();
}

double SimpleSpectrum::calculate_rate(double old_rate)
{
    return old_rate;
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
