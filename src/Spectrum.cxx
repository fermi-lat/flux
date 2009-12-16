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


double Spectrum::flux (double time ) const {
    return 0.; // flag that we don't have a flux
}

double Spectrum::solidAngle( )const
{
    return 1.0; // flag that doesn't calculate
}

std::pair<double,double> Spectrum::dir(double energy)
{
    // Purpose: return solid angle pair (costh, phi) for the given energy
    // Input:: the given energy.

#if 0 // disable this to allow Sun, Moon sources
    throw std::runtime_error("Spectrum::dir called: sub class must implement if use_sprectrum invoked");
#endif
    return std::make_pair(0,0);
}


const char * Spectrum::particleName()const{
  return m_particle_name.c_str();
}

double Spectrum::energy( double time)
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

double Spectrum::interval (double time)
{
    return -1.; //if this gets called, FluxSource will use the spectrum's flux.
                //to determine the interval.
}
