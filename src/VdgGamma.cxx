/** 
* @file VdgGamma.cxx
* @brief declaration and definition of VdgGamma
*
*/
#include "flux/Spectrum.h"
#include "flux/SpectrumFactory.h"

#include "CLHEP/Random/RandBreitWigner.h"
#include "CLHEP/Random/RandFlat.h"
#include <string>
#include <utility>
#include <algorithm>
#include <map>

/** 
* \class VdgGamma
*
* \brief Spectrum representing gammas produced by Vandegraaff accelerator at SLAC 
* \author Xin Chen
* 
*/
//

class VdgGamma : public Spectrum
{
public:

  VdgGamma(const std::string& paramstring) { }

  double energy( double time);
    
  virtual std::string title() const{return "VdgGamma";}
  virtual const char * particleName() const { return "gamma"; }
  inline  const char * nameOf() const {return "VdgGamma";}
    
    
private:

};


static SpectrumFactory<VdgGamma> factory;
const ISpectrumFactory& VdgGammaFactory = factory;

// VDG produces two typesof gammas:
// 17.619 MeV gammas with a width of 0.01 MeV, simulated as a line spectrum
// 14.586 MeV gammas with a width of 1.5 MeV, simulated as a BreitWigner spectrum 
// ratio between these two types of events is 100:49.7
double VdgGamma::energy( double time )
{
    static double ratio = 49.7/(49.7+100.);
    double r = RandFlat::shoot();

    if(r < ratio) {
      for(; ;) {
	float ene = RandBreitWigner::shoot(0.014586, 0.0015);

	//arbitrary 3 sigma cut
	static double low = 0.014586 - 3. * 0.0015;
	static double high = 0.014586 + 3. * 0.0015;
	if(ene > low && ene < high) return ene;
      }
    }
    else {
      return 0.0171619;
    }

}

  

