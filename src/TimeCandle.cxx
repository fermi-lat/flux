/**
 * @file TimeCandle.cxx
 * @brief Implementation of class TImeCandle.cxx: a source that ticks 

 * $Header$
 */

#include "TimeCandle.h"

#include "flux/FluxException.h" // for FATAL_MACRO
#include <utility>
#include <sstream>
#include <cmath>
#include "flux/SpectrumFactory.h"

static SpectrumFactory<TimeCandle> factory;
const ISpectrumFactory& TimeCandleFactory = factory;

TimeCandle::TimeCandle()
: m_T0(30.)
, m_name("TimeTick")
{}//default constructor
TimeCandle::TimeCandle(const std::string& params)
: m_T0(parseParamList(params,0)) 
, m_name("TimeTick")
, m_first(true)
{}



std::string TimeCandle::title()const
{
    std::stringstream s;
    s << "TimeTick("<< m_T0 << ")" ;
    std::string t(s.str()); 
    return t;
}



double TimeCandle::energy( double time)
{     
    m_first=false;  
    return 0.;
}


const char*
TimeCandle::particleName() const
{
    return m_name.c_str();
}

float TimeCandle::parseParamList(std::string input, int index)
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

double TimeCandle::interval (double time)
{  
    if( m_first){return 1e-30;} // epsilon to be greater than zero
    return m_T0;
}

