// $Id: AlbedoPSpectrum.cxx,v 1.11 2003/02/23 02:08:22 burnett Exp $


#include "AlbedoPSpectrum.h"

#include <cmath>
#include <algorithm>
#include <functional>
#include "CLHEP/Random/Random.h"
#include "flux/GPS.h"
#include "Geomag.h"
#include "flux/SpectrumFactory.h"

static SpectrumFactory<AlbedoPSpectrum> factory;
const ISpectrumFactory& AlbedoPSpectrumFactory = factory;


///Initializes parameters during construction
void AlbedoPSpectrum::init(const std::vector<float>& params) {
    // if there are parameters passed to chime, it shouuldn't move.
    m_allowMove = params.size()>1? false: true;
    float lat =  params.size()>0? params[0]: 0.0f;
    float lon =  params.size()>1? params[1]: 0.0f;
    setPosition(lat, lon);
    
    m_particle_name = "p";
    
    // set callback to be notified when the position changes
    m_observer.setAdapter( new ActionAdapter<AlbedoPSpectrum>
        (this, &AlbedoPSpectrum::askGPS) );
    
    GPS::instance()->notification().attach( &m_observer );
    
}


//-------------------------- constructors---------------------

///Constructor.  Initial geographic lat and lon in degrees.
AlbedoPSpectrum::AlbedoPSpectrum(const std::string& paramstring) {
    std::vector<float> params;
    
    parseParamList(paramstring,params);
    
    init(params);
}


//------------------------ title()

///What's my name?
std::string AlbedoPSpectrum::title() const {
    return "AlbedoPSpectrum";
}

//------------------------- particleName

///What kind of particle do I make?  (p for proton)
const char * AlbedoPSpectrum::particleName() const {
    return m_particle_name.c_str();
}

//-------------------------- setParticleName

///If anyone's crazy enough to override the default particle type.
void AlbedoPSpectrum::setParticleName(std::string name)
{
    m_particle_name = name;
}

//-------------------------- flux() (current position)

/// calculate flux for the current position
double AlbedoPSpectrum::flux(double) const {
    return m_flux;
}

//--------------------------

///The normalization used for the flux calculation.
double AlbedoPSpectrum::solidAngle() const
{
    return 4.*M_PI;
}

//-------------------------- flux() (specified position)

///Calculate flux for a specified position.
/// lat and lon are in degrees
///A simple formula that's consistent with the AMS data.
float AlbedoPSpectrum::flux(float lat, float lon) const {
    double v1,v2, alf1, alf2, emin, emax, Ejoin;
    fitParams(lat, lon, alf1, alf2, emin, emax, v1, v2, Ejoin);
    // The factor 1.47 is induced by the assumed angular dependence.
    return 1.47 * (v1+v2);
}

//-------------------------- flux() (position given in a pair)

///Flux for a specified position, packaged as a pair.
float AlbedoPSpectrum::flux( std::pair<double,double> coords) const {
    return flux(coords.first, coords.second);
}

//-------------------------- operator()  sample an energy value

/// Return a random value of energy sampled from the spectrum.
/// A simple analytical formula that's consistent with the spectrum
/// measured by AMS.
float AlbedoPSpectrum::operator() (float x) const{
    double v1,v2, alf1, alf2, emin, emax, Ejoin;
    fitParams(m_lat, m_lon, alf1, alf2, emin, emax, v1, v2, Ejoin);
    double split = v1/(v1+v2);
    if (x > split) {
        return pow((pow(Ejoin,1.-alf2)+((1.-x)/(1.-split))*(pow(emax,1.-alf2)
            -pow(Ejoin,1.-alf2))), 1./(1.-alf2));
    } else {
        return pow((pow(emin,1.-alf1)+(x/split)*(pow(Ejoin,1.-alf1)
            -pow(emin,1.-alf1))), 1./(1.-alf1));
    }
}

//-------------------------- askGPS()

///Ask the GPS where we are located.
int AlbedoPSpectrum::askGPS()
{
    if(m_allowMove)setPosition(GPS::instance()->lat(), GPS::instance()->lon());
    return 0; // can't be void in observer pattern
}

//-------------------------- setPosition(separate coordinates)
/// Do the initialization necessary when moving to a new position
void AlbedoPSpectrum::setPosition(float lat, float lon) {
    m_lat = lat;
    m_lon = lon;
    m_flux = flux(lat, lon);
}

//-------------------------- setPosition (coordinates packaged in a pair)
/// Do the initialization necessary when moving to a new position
void AlbedoPSpectrum::setPosition( std::pair<double,double> coords) {
    AlbedoPSpectrum::setPosition(coords.first, coords.second);
}

//------------------------- calculate_rate()
///Ask GPS where we are and find the flux at that point.
double AlbedoPSpectrum::calculate_rate(double)
{
    return flux(GPS::instance()->lat(), GPS::instance()->lon());
}

//-------------------------- dir()
///Choose random direction for a particle.
///Zenith angle distribution according to 1+0.6*sin(zenith)
std::pair<double,double> AlbedoPSpectrum::dir(double /*energy*/)
{
    // Random particle direction
    
    float coszenith, earthazi,dens,v;
    const int try_max = 1000;
    static int max_tried = 0;
    int trial = 0;
    earthazi = 2. * M_PI * HepRandom::getTheGenerator()->flat();
    do {
        coszenith = 2. * HepRandom::getTheGenerator()->flat() - 1.;
        dens = 1. + 0.6 * sqrt(1.-coszenith*coszenith);
        v = 1.6 * HepRandom::getTheGenerator()->flat();
    } while (v > dens && trial++ < try_max);
    
    max_tried = std::max(trial, max_tried);
    
    return std::make_pair<double,double>(coszenith, earthazi);
}

//------------------------- fitParams()
/** Evaluate the coefficients of the broken power law albedo proton model.
Input: lat  Magnetic latitude in degrees
Output: alf1,alf2   Slopes of the two power laws
emin,emax  Artificial upper and lower cutoff energies
v1,v2   Integrated fluxes of the two power laws
Ejoin   Energy at which the two power laws meet

  flux (protons / m^s sec ster MeV) = 
  a1 * (e/e1)^(-alf1)    if emin < e < Ejoin
  a2 * (e/e2)^(-alf2>    if Ejoin < e < emax
  when e is measured in GeV  (following the AMS convention)
*/
void AlbedoPSpectrum::fitParams(const double lat, const double lon, 
                                double& alf1, double& alf2, double& emin, double& emax, 
                                double& v1, double &v2, double& Ejoin) const
{
    double theta = abs(Geomag::geolat(lat,lon)) * M_PI / 180.;
    double e1 = .2;  // pivot point of lower power law
    double e2 = 1.;  // pivot point of higher power law
    emin = .01;      
    emax = 10.;
    // normalization values of the two power laws
    double a1 = 0.142 + theta * (-.4809 + theta * .5517);
    double a2 = .0499 + theta * (-.3239 + theta * (.8077 + theta * (-.8800 + 
        theta * .3607)));
    alf1 = .4913 + theta * (2.017 + theta * (-.6941 - 1.49 * theta));
    alf2 = 2.85 - .875 * theta;
    Ejoin = pow(a1*pow(e1,alf1) / (a2*pow(e2,alf2)), 1./(alf1-alf2));
    v1 = 1000.*a1*pow(e1,alf1)*(pow(Ejoin,1.-alf1)-pow(emin,1.-alf1))
        /(1.-alf1);  // integral of lower power law
    v2 = 1000.*a2*pow(e2,alf2)*(pow(emax,1.-alf2)-pow(Ejoin,1.-alf2)) 
        /(1.-alf2);  // integral of upper power law
}
