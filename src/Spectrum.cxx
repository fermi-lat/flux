/** @file Spectrum.cxx
    @brief  implementation of the Spectrum class.

    $Header: /nfs/slac/g/glast/ground/cvs/flux/src/Spectrum.cxx,v 1.8 2007/01/23 16:24:44 burnett Exp $
*/
#include "flux/Spectrum.h"
#include <cmath>
#include <cstdlib>
#include <stdexcept>

// CLHEP
#include "CLHEP/Random/RandFlat.h"

double Spectrum::s_startTime(0);


Spectrum::~Spectrum(){}


auto Spectrum::flux (double time ) const -> double {
    return 0.; // flag that we don't have a flux
}

auto Spectrum::solidAngle( )const -> double
{
    return 1.0; // flag that doesn't calculate
}

auto Spectrum::dir(double energy) -> std::pair<double,double>
{
    // Purpose: return solid angle pair (costh, phi) for the given energy
    // Input:: the given energy.

#if 0 // disable this to allow Sun, Moon sources
    throw std::runtime_error("Spectrum::dir called: sub class must implement if use_sprectrum invoked");
#endif
    return std::make_pair(0,0);
}


auto Spectrum::particleName()const -> const char *{
  return m_particle_name.c_str();
}

auto Spectrum::energy( double time) -> double
{
    // default implementation, which calls the operator()(float r)
    return (*this)(CLHEP::RandFlat::shoot());
}


void Spectrum::parseParamList(std::string input, std::vector<float>& output) const
{
    
    int i=0;
    for(;!input.empty() && i!=std::string::npos;){
        float f = ::atof( input.c_str() );
        output.push_back(f);
        i=input.find_first_of(",");
        input= input.substr(i+1);
    } 
}

auto Spectrum::interval (double time) -> double
{
    return -1.; //if this gets called, FluxSource will use the spectrum's flux.
                //to determine the interval.
}
