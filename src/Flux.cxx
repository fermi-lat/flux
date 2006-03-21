/** @file Flux.cxx
@brief Implementation of Flux

$Header: /nfs/slac/g/glast/ground/cvs/flux/src/Flux.cxx,v 1.8 2005/06/15 21:42:16 burnett Exp $

Original author: T. Burnett
*/
#include "flux/Flux.h"
#include "flux/EventSource.h"
#include "flux/FluxMgr.h"
#include "flux/SpectrumFactory.h"

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
std::string Flux::name() const
{
    return m_flux->name();
}

/// full title of the flux
std::string Flux::title()const 
{
    return m_event!=0? m_event->fullTitle() : "?";
}


void Flux::generate()
{
    // Purpose and Method: generate a new entry trajectory, set FluxSource, time
    // Inputs  - none
    // Outputs - none
    do{
        double current_time=time();
        // get the next event and its time interval?
        m_flux = m_event->event(current_time);
        double timepass = m_event->interval(current_time);
        pass(timepass);
    }while(m_event->occulted());
}

// the particle generated 
std::string Flux::particleName()const{
    return std::string(m_flux->particleName());
}

// its kinetic energy
double Flux::energy()const
{
    return m_flux->energy();
}

// starting point 
HepPoint3D Flux::launchPoint()const
{
    return m_flux->launchPoint();
}

double Flux::time()const 
{
    return s_mgr->time();
}


/// pass a specific amount of time    
void Flux::pass ( double t){
    s_mgr->pass(t);
}

/// Get the time as held by GPS    
double Flux::gpsTime () const{
    return s_mgr->time();
}


// direction
HepVector3D Flux::launchDir()const
{
    return m_flux->launchDir();
}

// rate ( /mm**2 /s)
double Flux::rate()const
{
    return  m_event!=0?  m_event->rate(time()) : -1;;
}

/// set the area of the target
void Flux::setTargetArea( double area)
{
    m_event->totalArea(area);
}

double Flux::targetArea()const
{
    return m_event->totalArea();
}


/// find which spectrum created the current particle
std::string Flux::findSource()const
{
    return m_event->findSource();
}

/// return a unique number correcponding to that spectrum
int Flux::numSource()const
{
    return m_event->numSource();

}


void Flux::addFactory(std::string name, const ISpectrumFactory* factory ) {
    SpectrumFactoryTable::instance()->addFactory(name,factory);
}

CLHEP::HepRotation Flux::CELTransform(double time)const{
    return s_mgr->CELTransform(time);
}
CLHEP::HepRotation Flux::orientTransform(double time)const{
    return s_mgr->orientTransform(time);
}

CLHEP::HepRotation Flux::transformGlastToGalactic(double time)const{

    return s_mgr->transformGlastToGalactic(time);
}

//get the transformtation matrix - the rest of these functions are now deprecated
CLHEP::HepRotation Flux::transformToGlast(double seconds,astro::GPS::CoordSystem index)const{
    return s_mgr->transformToGlast(seconds, index);
}



void Flux::writeSourceCharacteristic(std::ostream& out){
    m_event->writeSourceCharacteristic(out);
}
