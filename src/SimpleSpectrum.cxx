/** @file SimpleSpectrum.cxx
    @brief definition of SimpleSpectrum

   $Header: /nfs/slac/g/glast/ground/cvs/flux/src/SimpleSpectrum.cxx,v 1.13 2007/08/29 17:40:15 burnett Exp $
*/


#include "flux/SimpleSpectrum.h"
#include "flux/SpectrumFactory.h"
#include "flux/FluxException.h" // for FATAL_MACRO

#include <xercesc/dom/DOMElement.hpp>
#include "xmlBase/Dom.h"
#include "facilities/Util.h"

#include <cstdlib>
#include <utility>
#include <sstream>
#include <cmath>
#include <map>
#include <stdexcept>

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
    ///@class ParMap 
    ///@brief local analysis of keyword string
    class ParMap{
    public:
        ParMap(std::string paramString)
        {
            facilities::Util::keyValueTokenize(paramString,",",m_tokenMap);
        }
        double value(std::string name) const{
            std::map<std::string,std::string>::const_iterator it = m_tokenMap.find(name);
            if( it==m_tokenMap.end() ) {
                throw std::invalid_argument("SimpleSPectrum: keyword "+name+" not found");
            }
            return std::atof(it->second.c_str());
        }
    private:
        std::map<std::string,std::string> m_tokenMap;
    };

}


// implement simple broken power law
SimpleSpectrum::SimpleSpectrum(const std::string& paramString)
: m_name("gamma")
, m_E0(10)
, m_index(2.0)
, m_index2(2.0)
, m_ebreak(0)
, m_emax(200000)
{
    bool qemax(false); 
    ParMap parmap(paramString);
    try {
        m_E0 = parmap.value("emin");
    } catch (...) { }
    try {
        m_emax = parmap.value("emax");
        qemax=true;
    } catch (...) {}
    try {
        m_index = parmap.value("gamma");
    } catch (...) {}
    try {
        m_index2 = parmap.value("gamma2");
        m_ebreak = parmap.value("ebreak");
    } catch (...) {}

    setup_power_law();
}


SimpleSpectrum::SimpleSpectrum(const 
                               XERCES_CPP_NAMESPACE_QUALIFIER DOMElement* xelem,
                               bool useGeV)
: m_useGeV(useGeV)
{
    m_name = xmlBase::Dom::getAttribute(xelem, "name").c_str();
    
    const XERCES_CPP_NAMESPACE_QUALIFIER DOMElement* spectrum = 
      xmlBase::Dom::findFirstChildByName(xelem, "*");
    
    std::string tagName = xmlBase::Dom::getTagName(spectrum);
    if( tagName == "power_law" ){
        m_E0 =      xmlBase::Dom::getDoubleAttribute(spectrum, "emin");
        m_emax = xmlBase::Dom::getDoubleAttribute(spectrum, "emax");
        m_index =xmlBase::Dom::getDoubleAttribute(spectrum, "gamma"); 
        m_ebreak = xmlBase::Dom::getDoubleAttribute(spectrum, "ebreak");
        m_index2 =xmlBase::Dom::getDoubleAttribute(spectrum, "gamma2");

        setup_power_law();
    }
    else if(tagName=="energy") {
        // single energy: no interpolation
        m_emax =m_E0 =xmlBase::Dom::getDoubleAttribute(spectrum, "e");
    }
    else if( tagName == "exponential") {
        m_E0 = xmlBase::Dom::getDoubleAttribute(spectrum, "exponential");
        m_index = xmlBase::Dom::getDoubleAttribute(spectrum,"exponent");
        m_emax = 100.0;
        m_index = 0.0;
        FATAL_MACRO("exponential spectral component not implemented yet");
    }
    else {
        std::cerr << "Unknown name: " << m_name << std::endl;
        FATAL_MACRO("Unknown particle spectrum!");
    }
}

void SimpleSpectrum::setup_power_law()
{
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
SimpleSpectrum::operator()(float f)
{
    if( m_emax==m_E0)   return m_E0;
    
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
    // @todo: throw explicit exception
    if( index < output.size() ){
        return output.at(index);
    }else{
        return 0;
    }

}
