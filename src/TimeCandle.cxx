/**
 * @file TimeCandle.cxx
 * @brief Implementation of class TImeCandle.cxx: a source that ticks 

 * $Header: /nfs/slac/g/glast/ground/cvs/flux/src/TimeCandle.cxx,v 1.7 2006/11/15 22:30:09 burnett Exp $
 */

#include "TimeCandle.h"

#include "flux/FluxException.h" // for FATAL_MACRO
#include <cstdlib>
#include <utility>
#include <sstream>
#include <cmath>
#include "flux/SpectrumFactory.h"

static SpectrumFactory<TimeCandle> factory;
const ISpectrumFactory& TimeCandleFactory = factory;

TimeCandle::TimeCandle()
: m_period(30.)
, m_name("TimeTick")
{}//default constructor
TimeCandle::TimeCandle(const std::string& params)
: m_period(parseParamList(params,0, 1)) 
, m_offset(parseParamList(params,1,-1))
, m_name("TimeTick")
, m_first(true)
{}



std::string TimeCandle::title()const
{
    std::stringstream s;
    s << "TimeTick("<< m_period << ")" ;
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

float TimeCandle::parseParamList(std::string input, int index, float defaultValue)
{
    std::vector<float> output;
    int i=0;
    
    for(;!input.empty() && i!=std::string::npos;){
        float f = ::atof( input.c_str() );
        output.push_back(f);
        i=input.find_first_of(",");
        input= input.substr(i+1);
    } 
    if( index>output.size()-1 )return defaultValue;
    return output[index];
}

double TimeCandle::interval (double time)
{  
    if( m_first){
        if( m_offset<0) return 1e-30;// no offset; epsilon to be greater than zero

        // determing time until next offset from current time
        double tic = ::floor(time), next(tic+m_offset-time);
        return next>0? next : next+1.;

    } 
    return m_period;
}

