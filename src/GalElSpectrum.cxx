// $Header: /nfs/slac/g/glast/ground/cvs/flux/src/GalElSpectrum.cxx,v 1.4 2005/06/15 21:42:16 burnett Exp $
// Original author: P. L. Nolan, pln@egret1.Stanford.EDU
//
#ifdef __GNUG__
#pragma implementation
#endif

#include "GalElSpectrum.h"
#include "astro/GPS.h"
using astro::GPS;
#include <cmath>
#include "CLHEP/Random/Random.h"

#include "flux/SpectrumFactory.h"

static SpectrumFactory<GalElSpectrum> factory;
const ISpectrumFactory& GalElSpectrumFactory = factory;


const float GalElSpectrum::m_rearth = 6371.f;
const float GalElSpectrum::m_altitude = 600.f;

GalElSpectrum::GalElSpectrum(const std::string& paramstring):m_pspec(paramstring) {
    std::vector<float> params;
    
    parseParamList(paramstring,params);

    // if there are parameters passed to chime, it shouldn't move.
    m_allowMove = params.size()>1? false: true;
    
    float lat = params.size()>0? params[0] : 0;
    float lon = params.size()>1? params[1] : 0;
    init(lat, lon); 
}


GalElSpectrum::~GalElSpectrum() {};

double GalElSpectrum::calculate_rate(double old_rate) {
    return flux(GPS::instance()->lat(), GPS::instance()->lon());
}

double GalElSpectrum::flux(double) const {
    return -m_normfact * (m_norm / m_expo) * pow(m_cutoff, m_expo);
}

double GalElSpectrum::solidAngle() const {
    return 4. * M_PI;
}

float GalElSpectrum::flux(float cut) const {
    return -m_normfact * (m_norm / m_expo) * pow(cut, m_expo);
}

float GalElSpectrum::flux(float lat, float lon) const {
    float Ecut = findCutoff(lat,lon);
    return flux(Ecut);
}

float GalElSpectrum::flux(std::pair<double,double> coords) const {
    return flux(coords.first, coords.second);
}

float GalElSpectrum::operator() (float x) {
    return m_cutoff * pow(1.-x, 1./m_expo);
}

void GalElSpectrum::setPosition(float lat, float lon) {
    m_lat = lat;
    m_lon = lon>0? lon : lon+360.;
    m_pspec.setPosition(m_lat, m_lon);
    m_cutoff = findCutoff(m_lat,m_lon);
    m_flux = flux(m_cutoff);
}

void GalElSpectrum::setPosition(std::pair<double,double> coords) {
    GalElSpectrum::setPosition(coords.first, coords.second);
}

int GalElSpectrum::askGPS() {
    if(m_allowMove)setPosition(GPS::instance()->lat(), GPS::instance()->lon());
    return 0;
}

float GalElSpectrum::findCutoff(float rflux) const {
    return pow(-m_expo*rflux/m_norm, (float) 1./m_expo);
}
float GalElSpectrum::findCutoff(float lat, float lon) const {
#if 0
    float ProtonCut = m_pspec.findCutoff(lat,lon);
#else
    float ProtonCut = m_pspec.findCutoff(); // from current
#endif
    return sqrt(ProtonCut*(ProtonCut+2.*0.938));
}
float GalElSpectrum::findCutoff(std::pair<double,double> coords) const {
    return findCutoff(coords.first, coords.second);
}

std::pair<double,double> GalElSpectrum::dir(double energy)  {
    float earthazi = 2.*M_PI* CLHEP::HepRandom::getTheGenerator()->flat();
    float coszenith = m_coscutoff + (1.-m_coscutoff)*
        CLHEP::HepRandom::getTheGenerator()->flat();
    return std::make_pair<float,float>(coszenith, earthazi);
}

std::string GalElSpectrum::title() const {
    return "GalElSpectrum";
}

const char * GalElSpectrum::particleName() const {
    return m_particle_name.c_str();
}

void GalElSpectrum::setParticleName(std::string name) {
    m_particle_name = name;
}

void GalElSpectrum::init(float lat, float lon) {
    m_normfact = .5*(1.+sqrt(m_altitude*(m_altitude+2.*m_rearth)) /
        (m_rearth+m_altitude));
    m_expo = static_cast<float>(-3.086 + 1.);
    m_norm = 227.f;
    setPosition(lat, lon);
    m_coscutoff = -sqrt(m_altitude*m_altitude+2.*m_altitude*m_rearth)
        / (m_altitude+m_rearth);
    m_particle_name = "e-";
    m_observer.setAdapter( new ActionAdapter<GalElSpectrum>(this,
					   &GalElSpectrum::askGPS) );
    GPS::instance()->notification().attach( &m_observer );
}




