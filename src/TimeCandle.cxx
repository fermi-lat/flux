// $Header: /nfs/slac/g/glast/ground/cvs/flux/src/TimeCandle.cxx,v 1.1.1.1 2003/07/29 18:22:19 burnett Exp $


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
{}



std::string TimeCandle::title()const
{
    std::stringstream s;
    s << particleName() << '(' << 1 << " GeV";
    s << ")" << '\0';
    std::string t(s.str()); 
    return t;
}

float
TimeCandle::operator()(float f)const
{
    return 1.;
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
