// Spectrum.cxx: implementation of the Spectrum class.
//
//////////////////////////////////////////////////////////////////////

#include "flux/Spectrum.h"
#include <cmath>

// CLHEP
#include "CLHEP/Random/RandomEngine.h"
#include "CLHEP/Random/RandGeneral.h"
#include "CLHEP/Random/RandExponential.h"
#include "CLHEP/Random/RandFlat.h"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

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
    // default: random except for Earth occultation
    //here's an attempt at random default distributions as above:
    return std::make_pair(((RandFlat::shoot(1.0))*1.4)-0.4,(RandFlat::shoot(1.0))*2*M_PI);
    
}


const char * Spectrum::particleName()const{
  return m_particle_name.c_str();
}

double Spectrum::energy( double time)
{
    // default implementation, which works for other Spectrum objects
    return (*this)(RandFlat::shoot());
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
