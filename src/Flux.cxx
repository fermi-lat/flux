/** @file Flux.cxx
@brief Implementation of Flux

$Header: /nfs/slac/g/glast/ground/cvs/flux/src/Flux.cxx,v 1.18 2008/01/07 04:18:23 burnett Exp $

Original author: T. Burnett
*/
#include "flux/Flux.h"
#include "flux/EventSource.h"
#include "flux/FluxMgr.h"
#include "flux/SpectrumFactory.h"
#include <stdexcept>

Flux::Flux(std::string name) 
:  m_flux(0)
{
    m_event = s_mgr->source(name);
}
Flux::Flux(std::vector<std::string> names) 
:  m_flux(0)
{
    m_event = s_mgr->compositeSource(names);
}
Flux::~Flux() 
{
    delete m_flux;
}

FluxMgr* Flux::s_mgr=0;

void Flux::mgr(FluxMgr* m){ s_mgr=m;}
// name of the flux
auto Flux::name() const -> std::string
{
    return m_flux->name();
}

/// full title of the flux
auto Flux::title()const -> std::string 
{
    return m_event!=0? m_event->fullTitle() : "?";
}


auto Flux::generate() -> bool
{
    // Purpose and Method: generate a new entry trajectory, set FluxSource, time
    // Inputs  - none
    // Outputs - valididty status
    do{
        double current_time=time();
        // get the next event and its time interval?
        if( ! m_event->enabled()) return false; // there is no source
        m_flux = m_event->event(current_time);
        double delta_time = m_event->interval();
        if( delta_time < 0) {
            throw std::runtime_error("Flux::generate:  negative interval generated");
        }
        setTime( current_time + delta_time);
        s_mgr->synch();  // notify observers
    }while(m_event->occulted() && m_flux->enabled());
    return m_flux->enabled();
}

// the particle generated 
auto Flux::particleName()const -> std::string{
    return std::string(m_flux->particleName());
}

// its kinetic energy
auto Flux::energy()const -> double
{
    return m_flux->energy();
}

// starting point 
auto Flux::launchPoint()const -> CLHEP::Hep3Vector
{
    return m_flux->launchPoint();
}

auto Flux::time()const -> double
{
    return s_mgr->time();
}

void Flux::setTime(double newtime)
{
    s_mgr->setTime(newtime);
}

#if 0
/// pass a specific amount of time    
void Flux::pass ( double t){
    s_mgr->pass(t);
}
/// Get the time as held by GPS    
double Flux::gpsTime () const{
    return s_mgr->time();
}
#endif


// direction
auto Flux::launchDir()const -> CLHEP::Hep3Vector
{
    return m_flux->launchDir();
}

// rate ( /mm**2 /s)
auto Flux::rate()const -> double
{
    return  m_event!=0?  m_event->rate(time()) : -1;;
}

/// set the area of the target
void Flux::setTargetArea( double area)
{
    m_event->totalArea(area);
}

auto Flux::targetArea()const -> double
{
    return m_event->totalArea();
}


/// find which spectrum created the current particle
auto Flux::findSource()const -> std::string
{
    return m_event->findSource();
}

/// return a unique number correcponding to that spectrum
auto Flux::numSource()const -> int
{
    return m_event->numSource();

}


void Flux::addFactory(std::string name, const ISpectrumFactory* factory ) {
    SpectrumFactoryTable::instance()->addFactory(name,factory);
}
//get the transformtation matrix - the rest of these functions are now deprecated
auto Flux::transformToGlast(double seconds,astro::GPS::CoordSystem index)const -> CLHEP::HepRotation{
    return s_mgr->transformToGlast(seconds, index);
}



void Flux::writeSourceCharacteristic(std::ostream& out){
    m_event->writeSourceCharacteristic(out);
}

auto Flux::invalid()const -> bool{
    return !m_flux->enabled();
}

