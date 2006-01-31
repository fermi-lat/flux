/** 
* @file SurfaceMuons.cxx
* @brief declaration and definition of SurfaceMuons
*
*  $Header: /nfs/slac/g/glast/ground/cvs/flux/src/SurfaceMuons.cxx,v 1.12 2005/03/27 03:01:46 burnett Exp $
*/
#include "flux/Spectrum.h"
#include "flux/SpectrumFactory.h"

#include "CLHEP/Random/RandFlat.h"
#include <string>
#include <utility>
#include <algorithm>
#include <map>
#include <cmath>
/** 
* \class SurfaceMuons
*
* \brief Spectrum representing cosmic ray muon flux at the Earth's surface
* \author T. Burnett
* 
* $Header: /nfs/slac/g/glast/ground/cvs/flux/src/SurfaceMuons.cxx,v 1.12 2005/03/27 03:01:46 burnett Exp $
*/
//

class SurfaceMuons : public Spectrum
{
public:
    /** @brief ctot
    @param paramstring string from xml. If present, assume is the range of costheta to generate
    */
    SurfaceMuons(const std::string& paramstring);

    /// flux integrated over energy, will be multipled by solidAngle to get a rate/m^2
    double flux (double /*time*/ ) const { return m_flux/solidAngle();}

    double solidAngle()const{return 2*M_PI*fabs(m_cosmax-m_cosmin);}

    /** @brief sample a single particle energy from the spectrum
    @param current time (ignored)
    @return the energy in GeV, or MeV if the source was specified as such
    */
    double energy( double time);


    /** @brief return solid angle pair (costh, phi) for the given energy
    @param energy the energy previously generated
    */
    virtual std::pair<double,double> dir(double energy);
    
    virtual std::string title() const{return "SurfaceMuons";}

    virtual const char * particleName() const;
    inline  const char * nameOf() const {return "SurfaceMuons";}


private:

  // option to choose which energy spectrum is implemented.
  // 0 (default): use the spectrum as a function of E*cos(theta)
  // 1: use Hiro's analytic form, which is independent of angle, see analyticSpectrum()
  // 2: use data from Caprice94 experiment (PRL 83, 21:  22 Nov. 1999)
  // 3: use data from Caprice97 experiment (PRL 83, 21:  22 Nov. 1999)
  int m_option;



    double analyticSpectrum(double time) const;

    double capriceSpectrum(double time) const;

    // local function that approximates the spectrum as a function of E*cos(theta)

    double spectrum(double e);

    /**
    This is the integral spectrum as a map where the key is the integral spectrum and the
    value the associated energy. This makes is easy to invert it by using the STL algorithm lower_bound.
    */
    std::map<double,double> m_ispec; // integral spectrum, for sampling
    double m_flux;
    double m_index;
    double m_emax;
    double m_cosmin, m_cosmax;
    double m_costh;
    double m_total;

};


static SpectrumFactory<SurfaceMuons> factory;
const ISpectrumFactory& SurfaceMuonsFactory = factory;

SurfaceMuons::SurfaceMuons(const std::string& paramstring)
: m_index(2.71)
, m_emax(1000)
, m_total(0)
{
    //Purpose: Initializes parameters during construction

    std::vector<float> params;

    parseParamList(paramstring,params);

    // determine integral flux for specified costheta range

    m_cosmin = (params.size()>0)? params[0]: 0.0;
    m_cosmax = (params.size()>1)? params[1]: 1.0;

    // total flux is (above 1 GeV) the integral over the cos theta, assuming cos(theta)**2
    // and PDG value of 70/m^2/s/sr in the vertical
    m_flux = 70. * 2*M_PI * fabs(pow(m_cosmax,3) - pow(m_cosmin,3))/3;

    // option > 0  is for alternate version
    m_option = (int) (params.size()>2 ? params[2]: 0.);

    if(m_option == 0) {

        // create integral table of the flux function, as a map of
        // energy and e*flux(e), with logarithmic energies
        int n=100;

        for( int i=0; i< n; ++i){
            double 
                e = pow(10., 0.025*i),
                f= e*spectrum(e);

            m_ispec[m_total +=f] =e;
        }
    } else {
        m_flux *= 100./70.;
    }
}


const char* SurfaceMuons::particleName()const
{
    /// purpose: return a point to the particle name, either mu+ or mu-
    static const char * pnames[] = {"mu+", "mu-"};
    static double 
        charge_ratio = 1.2,  // from many measurements.
        plus_fraction=1/(1+charge_ratio);
    return CLHEP::RandFlat::shoot()>plus_fraction? pnames[0]:pnames[1];
}
namespace {
    inline double cube(double x){return x*x*x;}
}

/// sample a single particle energy from the spectrum: assume called first
double SurfaceMuons::energy( double time )
{
    using namespace std;
    // first choose the angle for the dir function, assuming cos**2 distribution
    static double third=1.0/3.0;
    double r = CLHEP::RandFlat::shoot();
    m_costh = pow(r*cube(m_cosmax)+(1-r)*cube(m_cosmin), third);

    if(m_option == 1) return analyticSpectrum(time);
    if(m_option == 2 || m_option == 3) return capriceSpectrum(time);


    // select an energy by inverting the integral distribution
    double energy = 0; 

    double partial = CLHEP::RandFlat::shoot()*m_total; // the value of the integral to find the energy
    map<double,double>::const_iterator 
        element  = m_ispec.lower_bound(partial);

    if( element == m_ispec.begin()){
        // this should not happen, but catch it anyway
        energy = element->second;
    }	

    else if( element!=m_ispec.end() ){
        // interpolate here
        map<double,double>::const_iterator prev = element;
        prev--;
        double 
            dI = partial - prev->first, // the part of the integral
            deltaI = element->first - prev->first, //total increment
            energy_prev = prev->second, 
            deltaE = element->second - energy_prev; // totlal energy diff

        energy =  energy_prev + dI*deltaE/deltaI;
    } else {
        // at the top
        energy = m_emax;
    }

    // scale by cos(theta)
    return energy/m_costh;
}


double SurfaceMuons::capriceSpectrum(double time) const
{
   // Data taken from J. Kremer, et. al., Measurements of Ground-Level Muons at Two Geomagnetic Locations
   // PRL Vol 83, Num. 21, 22 November 1999

   // Kinetic energy in GeV
   double ke[] = {1.2054e-001,  2.1240e-001,  3.0806e-001,  4.5440e-001,  6.0227e-001,  7.5088e-001,  8.9991e-001,
                  1.0990e+000,  1.2983e+000,  1.4978e+000,  1.9970e+000,  2.8362e+000,  4.0157e+000,  5.3954e+000,  
                  6.8951e+000,  9.8949e+000,  1.5395e+001,  2.2895e+001,  3.0995e+001,  4.3494e+001,  6.0994e+001,
                  8.5494e+001,  1.1989e+002};


   // Integrated mu+ and mu- flux from Caprice94 data (Manitoba, Canada at atmospheric depth of 1000 g/cm^2)
   // (m^2 sr s)^-1
   double integ_flux94[] = {2.5000e+000,  5.5400e+000,  1.0280e+001,  1.4795e+001,  1.9130e+001,  2.3165e+001, 
      2.7965e+001,  3.2385e+001,  3.6165e+001,  4.4365e+001,  5.4277e+001,  6.3127e+001,
      6.9351e+001,  7.3641e+001,  7.8411e+001,  8.2250e+001,  8.4095e+001, 8.4978e+001,
      8.5565e+001,  8.5863e+001,  8.6024e+001,  8.6114e+001};

   // Integrated mu+ and mu- flux from Caprice97 data (New Mexico at atmospheric depth of 886 g/cm^2)
   // Note:  One measurement was missing in the 10th bin for mu+ in the paper (corresponding to p=1.84 GeV/c.  
   // The mu- value was multiplied by 1.2 to obtain a guess for the mu+ value. This value
   // was then used used in generating this array.
   double integ_flux97[] = {2.2700e+000,  5.4900e+000,  1.1220e+001,  1.7190e+001,  2.2890e+001,  2.8230e+001,
      3.4730e+001,  4.0390e+001,  4.5450e+001,  5.5550e+001,  6.6218e+001,  7.6437e+001,
      8.3682e+001,  8.8407e+001,  9.3807e+001,  9.7844e+001,  9.9869e+001,  1.0075e+002,
      1.0136e+002,  1.0167e+002,  1.0185e+002,  1.0194e+002 };

   double *integ_flux = (m_option == 2)?  integ_flux94 : integ_flux97;

   double target = CLHEP::RandFlat::shoot(integ_flux[21]);
   double m = 0, b = 0;
   
   for(unsigned int i = 0; i < 22; i++) {
      if(integ_flux[i] >= target) {
         // interpolate between bins
         if(i != 0) 
            m = ( ke[i+1] - ke[i] ) / ( integ_flux[i] - integ_flux[i-1] );              
         else 
            m = ( ke[i+1] - ke[i] ) / integ_flux[i];

         b = ke[i+1] - m * integ_flux[i];
         break;
      }
   }
    
   return (m * target + b)/m_costh;

}


std::pair<double,double> SurfaceMuons::dir(double)
{
    // purpose: return the pair; note uses the previous value for cos(theta)
    return std::make_pair(m_costh, CLHEP::RandFlat::shoot(2*M_PI));
}

static inline double sqr(double x){return x*x;}

double SurfaceMuons::spectrum(double ecth)
{
    // purpose: spectrum as a function of E*cos(theta)
    // returns value differential in E, according to stuff in particle properties

    double 
        atmos = 1/(1+1.1*ecth/115.) + 0.054/(1+1.1*ecth/850.),
        cutoff = ecth<30? exp(-sqr(::log10(ecth/30.))/0.55) : 1.0;


    return ecth<1.0 ? 0 : pow(ecth, -2.71)*atmos*cutoff;
}

double SurfaceMuons::analyticSpectrum(double time) const 
{
    // analytic form is based on following energy spectrum
    //  0.003 flat (0.2 GeV < E < 1. GeV)
    //  0.003 * E^(-1.1) (1 GeV < E < 4 GeV)
    //  0.006 * E^(-1.6) (4 GeV < E < 10 GeV)
    //  0.038 * E^(-2.6) (10 GeV < E < 200 GeV)
    // only generated events with energy from 0.2 GeV to 200 GeV

    // ratio used later in determining in which energy range the event is going
    // to be generated
    double norm = 0.0024 + 0.00388 + 0.00184 + 0.000592;
    double ratio[] = {0.0024/norm, 0.00388/norm, 0.00184/norm, 0.000592/norm};
    double energy;

    static CLHEP::HepRandom* randomEng = CLHEP::HepRandom::getTheGenerator();

    // determine the energy range of the event to be generated
    float range = randomEng->flat();

    // use transformation method to generate random numbers
    float rand = randomEng->flat();

    if(range <= ratio[0]) {  // 0.2 GeV < E < 1 GeV
        energy = 0.2 + 0.8 * rand;
    }
    else if(range <= ratio[0] + ratio[1]) {  // 1 GeV < E < 4 GeV
        static double normE0 = 1. - pow(4., -0.1);
        energy = pow(1. - normE0 * rand, -10.);
    }
    else if(range <= ratio[0] + ratio[1] + ratio[2]) {  // 4 GeV < E < 10 GeV
        static double normE1 = pow(4., -0.6) - pow(10., -0.6);
        static double temp = pow(4., -0.6);
        energy = pow(temp - normE1 * rand, -5./3.);
    }
    else { // 10 GeV < E < 200 GeV
        static double normE2 = pow(10., -1.6) - pow(200., -1.6);
        static double temp = pow(10., -1.6);
        energy = pow(temp - normE2 * rand, -1./1.6);
    }

    return energy;
}
