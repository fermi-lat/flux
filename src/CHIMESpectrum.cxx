// $Header: /nfs/slac/g/glast/ground/cvs/flux/src/CHIMESpectrum.cxx,v 1.1.1.1 2003/07/29 18:22:19 burnett Exp $


#include "CHIMESpectrum.h"

#include "CHIMESpectrum.inc"  // numerical data: energies, fluxes, gfluxes

#include <cmath>
#include <algorithm>
#include <functional>
#include "CLHEP/Random/Random.h"
#include "flux/GPS.h"

// this is needed to include in the executable or dll
#include "flux/SpectrumFactory.h"

static SpectrumFactory<CHIMESpectrum> factory;
const ISpectrumFactory& CHIMESpectrumFactory = factory;


// First deal with subclass CHIMESpectrum::InterpVec::
//  which is a vector<float> with a couple of methods added for
//  searching and interpolation.
//-------------------------- InterpVec constructor

CHIMESpectrum::InterpVec::InterpVec() : std::vector<float>() {}
// default constructor.

//-------------------------- InterpVec::search

CHIMESpectrum::Intrp CHIMESpectrum::InterpVec::search(float x) const {
    // Purpose: Binary search a vector for values that straddle x using STL lower_bound.
    // Deal gracefully with values off the end by extrapolation.
    // Contents of vector must be monotonic.
    // Return integer index for use by another vector
    
    std::vector<float>::const_iterator loc;  // search. ascending or descending?
    if (back() > front()) loc=std::lower_bound(begin(), end(), x, std::less<float>());
    else                  loc=std::lower_bound(begin(), end(), x, std::greater<float>());
    // lower_bound points to the entry after x
    
    int ind;  // convert STL iterator to integer index
    if (loc == begin())    ind = 1;   // extrapolate below table
    else if (loc == end()) ind = (end() - begin()) - 1;  // above table
    else                   ind = (loc - begin());  // in the table
    
    return Intrp(ind, (x-(*this)[ind-1]) / ((*this)[ind]-(*this)[ind-1]));
}

//-------------------------- InterpVec::interpolate

float CHIMESpectrum::InterpVec::interpolate(Intrp y) const {
    // linear interpolation between two vector values
    return (*this)[y.first-1] +
        y.second * ((*this)[y.first]-(*this)[y.first-1]);
}

//-------------------------- Beginning of CHIMESpectrum proper

const float CHIMESpectrum::m_rearth = 6371.f;  // radius of earth in km
const float CHIMESpectrum::m_altitude = 600.f;  // altitude of circular orbit


void CHIMESpectrum::init(std::string paramstring) {
    //Purpose: Initializes parameters during construction
    
    std::vector<float> params;
    
    parseParamList(paramstring,params);
    
    int nen = sizeof(energies)/sizeof(float);
    m_en.reserve(nen);
    int i;
    for (i=0; i< nen; i++) m_en.push_back(energies[i]/1000.);
    // convert MeV to GeV
    m_normfact = .5*(1.+sqrt(m_altitude*(m_altitude+2.*m_rearth)) /
        (m_rearth+m_altitude));
    //geometrical shadow factor in 600 km orbit
    m_fluxes.reserve(nen);
    m_fl.reserve(nen);
    
    const float width = 0.115f; // CHIME log spacing between standard energies
    m_etop = m_en.back() * (1.+.5*width);  // boundary between table & tail
    m_expo = -2.75 + 1.0;  // power law exponent (integral) of high-energy tail
    
    // Populate table of differential free space fluxes -- no geomagnetic cutoff
    float tfl = 0.f;
    for (i=74; i >= 0; i--) {
        tfl += width*1000.*m_en[i]*m_normfact*fluxes[i];
        m_fluxes.insert(m_fluxes.begin(), tfl);
    }
    
    // table of total flux as a function of latitude and longitude
    for (int ii=0; ii < 73; ii++) {
        for (int jj=0; jj < 13; jj++) {
            m_fluxTbl[ii][jj] = gfluxes[jj+13*ii];
        }
    }
    // if there are parameters passed to chime, it shouldn't move.
    m_allowMove = params.size()>1? false: true;
    // set the initial location
    double lat =  params.size()>0? params[0]: 0.0f;
    double lon =  params.size()>1? params[1]: 0.0f;
    
    setPosition(lat, lon);
    
    // cos(angle between zenith and horizon)
    m_coscutoff = -sqrt(m_altitude*m_altitude+2.*m_altitude*m_rearth)
        / (m_altitude+m_rearth);
    
    m_particle_name = "p";
    
    // set callback to be notified when the position changes
    m_observer.setAdapter( new ActionAdapter<CHIMESpectrum>(this,&CHIMESpectrum::askGPS) );
    
    GPS::instance()->notification().attach( &m_observer );
    
}


//-------------------------- constructors---------------------

CHIMESpectrum::CHIMESpectrum(const std::string& params) {
    init(params);
}


std::string CHIMESpectrum::title() const {
    return "CHIMESpectrum";
}

const char * CHIMESpectrum::particleName() const {
    return m_particle_name.c_str();
}
void CHIMESpectrum::setParticleName(std::string name)
{
    m_particle_name = name;
}

//-------------------------- flux()  (specified cutoff value)

float CHIMESpectrum::flux(float cut) const {
    // Total flux in protons / m^2 sec ster
    // Interpolate in table if possible, otherwise assume power law
    //  tail at high energy.
    if (cut > m_etop) return m_upper * pow(cut/m_etop, m_expo);
    else              return m_upper + m_fl.interpolate(m_en.search(cut));
}

//-------------------------- flux() (current cutoff value)

double CHIMESpectrum::flux(double) const {
    // calculate flux for the current cutoff
    return m_flux;
}

//--------------------------

double CHIMESpectrum::solidAngle() const
{
    // the normalization for the flux calculation (according to doc)
    return 4.* M_PI;
}

//-------------------------- flux() (specified position)

float CHIMESpectrum::flux() const {
    // Flux as a function of latitude and longitude in a 600 km orbit.
    // Linear interpolate in a table with a 5 degree sampling grid.
    double lat = m_lat, lon=m_lon;
    int ilat = static_cast<int>(lat/5.+6.);
    double/*float*/ a = fmod(lat+30., 5.)/5.;
    int ilon = static_cast<int>(lon/5.);
    double/*float*/ b = fmod(lon, 5.)/5.;
    return m_fluxTbl[ilon][ilat] * (1.-a) * (1.-b) +
        m_fluxTbl[ilon+1][ilat] * (1.-a) * b +
        m_fluxTbl[ilon][ilat+1] * a * (1.-b) +
        m_fluxTbl[ilon+1][ilat+1] * a * b;
}


//-------------------------- operator()  sample an energy value

float CHIMESpectrum::operator() (float x)const {
    // return a random value of energy sampled from the spectrum
    
    float join = (m_tot-m_upper)/m_tot;
    if (x < join) return m_en.interpolate(m_fl.search((1.-x)*m_tot-m_upper));
    else          return m_etop*pow((1.-x)/(1.-join), 1./m_expo);
}

//-------------------------- setPosition (separate coordinates)

int CHIMESpectrum::askGPS()
{
    if(m_allowMove)setPosition(GPS::instance()->lat(), GPS::instance()->lon());
    return 0; // can't be void in observer pattern
}

void CHIMESpectrum::setPosition(double lat, double lon) {
    // Do the initialization necessary when moving to a new position:
    // look up cutoff energy, build a new table of integral proton
    // fluxes

    m_lat = lat;
    m_lon = lon>0? lon : lon+360.;
    
    // Integrated flux in the power law tail above the table.
    m_upper = -0.115*1000.*m_en.back()*m_normfact*fluxes[74]
        * pow((float)(1.+.5*0.115),m_expo)/(0.115*m_expo);
    
    m_cutoff = findCutoff();
    
    // Populate table of integral fluxes modified by geomagnetic cutoff.
    float tfl = 0.;
    m_fl.erase(m_fl.begin(), m_fl.end());
    for (int i = 74; i >= 0; i--) {
        tfl += 0.115*1000.*m_en[i]*m_normfact*fluxes[i]*exposure(m_en[i]);
        m_fl.insert(m_fl.begin(), tfl);
    }
    
    m_tot = m_fl.front() + m_upper;
    m_flux = flux(m_cutoff);
}

//-------------------------- setPosition (coordinates packaged in a pair)

void CHIMESpectrum::setPosition(std::pair<double,double> coords) {
    CHIMESpectrum::setPosition(coords.first, coords.second);
}
//-------------------------- findCutoff (from the current position)

float CHIMESpectrum::findCutoff() const {
    // determine the cutoff value at the geographical location
    return findCutoff(flux());
}

//-------------------------- findCutoff (from a flux)

float CHIMESpectrum::findCutoff(float rflux) const {
    // determine the cutoff value which will produce the desired flux
    return m_en.interpolate(m_fluxes.search(rflux-m_upper));
}


//------------------------- rad2()

float CHIMESpectrum::rad2() const {
    // square of (distance from center of magnetic dipole / earth radius)
    // Dipole is offset from the earth's center
    
    double/*float*/ d2 =
        pow((m_rearth+m_altitude)*sin(m_lat)            - 145.1, 2) +
        pow((m_rearth+m_altitude)*cos(m_lat)*cos(m_lon) + 371.2, 2) +
        pow((m_rearth+m_altitude)*cos(m_lat)*sin(m_lon) - 233.7, 2);
    return d2 / (m_rearth*m_rearth);
}

//------------------------- cosomega()

float CHIMESpectrum::cosomega(float E) const {
    // Opening angle of geomagnetic cutoff cone.
    // This is a pretty backward way to do it.  It starts from a cutoff
    // energy which is derived from a desired rate.  Then it calculates
    // the magnetic latitude from that.  Really, the cutoff and latitude
    // should be derived from the position by looking in tables.
    // Also, this is simple Størmer theory, ignoring  the penumbral
    // region and multipole effects.

    E *= 0.001;  // Convert Mev->GeV
    
    const float Mp = 0.938f;  // mass of proton in GeV
    double/*float*/ pcut = sqrt(m_cutoff*m_cutoff + 2.*Mp*m_cutoff);  // cutoff momentum
    const float moment = 59.8f; // magnetic moment of earth in convenient units
    double/*float*/ coslambda4 = 4. * pcut  * rad2() / moment;  // magnetic latitude**4
    double/*float*/ p = sqrt(E*E + 2.*Mp*E); // momentum
    double/*float*/ coso = 4. * (sqrt(pcut/p) - pcut/p) * pow(coslambda4, -0.75);
    // opening angle of Størmer cone
    if (coso > 1.) return 1.;
    else if (coso < -1.) return -1.;
    else return coso;
}

//-------------------------- exposure()

float CHIMESpectrum::exposure(float E) const {
    // Geomagnetic cutoff fraction.  Varies from 0 to 1.
    return 0.5 * (1. + cosomega(E));
}

//-------------------------- dir()

std::pair<double,double> CHIMESpectrum::dir(double energy)
{
    
    // Random particle direction from Stormer cone
    
    // Rejection method for the direction.  Direction must be inside
    // the Stormer cone and above the horizon.
    double/*float*/ coszenith, sinpolar, cospolar, azi;
    const int try_max = 1000;
    static int max_tried = 0;
    int trial=0;
    do {
        // uniform distribution within Stormer cone
        cospolar = -1. + HepRandom::getTheGenerator()->flat()*(cosomega(energy)+1.);
        sinpolar = sqrt(1.-cospolar*cospolar);
        azi = 2.*M_PI* HepRandom::getTheGenerator()->flat();
        coszenith = cos(azi)*sinpolar;
    } while (coszenith < m_coscutoff && trial++ < try_max);
    
    max_tried = std::max(trial, max_tried); // keep track
    
    // Transform to local earth-based coordinates.
    float earthazi = atan2(sinpolar*sin(azi), cospolar);
    
    return std::make_pair<double,double>(coszenith, earthazi);
}
