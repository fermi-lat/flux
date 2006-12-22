/** @file Flux.h
    @brief Declaration of Flux

   $Header: /nfs/slac/g/glast/ground/cvs/flux/flux/Flux.h,v 1.9 2006/11/07 03:34:28 burnett Exp $

  Original author: T. Burnett
*/

#ifndef FLUXSVC_FLUX_H
#define FLUXSVC_FLUX_H

/** 
* \class Flux
*
* \brief The class holding the interface with FluxMgr, EventSource, and FluxSource of the flux package. 
* Flux is used to get the actual information(energy, name, etc) about the current particle, and to generate new ones,
* through this interface.
* \author Toby Burnett
* 
* $Header $
*/

#include "IFlux.h"
#include "CLHEP/Vector/Rotation.h"

#include <vector>

// forward declarations
class FluxMgr;
class EventSource;

class Flux : public IFlux {
public:
    /// ctor, select the name
    Flux(std::string name);

    Flux(std::vector<std::string> names);

    virtual ~Flux();
    
    /// name of the flux
    virtual std::string name()const;
    
    /// full title of the flux
    virtual std::string title()const;
    
    /// generate a new entry trajectory
    virtual bool generate();
    
    /// the particle generated 
    virtual std::string particleName()const;
    
    /// the particle property entry for the last particle generated 
    //virtual ParticleProperty* property()const=0;
    
    /// its kinetic energy
    virtual double energy()const;
    
    /// starting point 
    virtual Hep3Vector launchPoint()const;
    
    /// direction
    virtual Hep3Vector launchDir()const;
    
    /// return the time
    virtual double time()const;
    
    /// pass a specific amount of time
    virtual void pass ( double t);
    
    /// Get the time as held by GPS (same: here for backward compatibility
    double gpsTime () const;
    
    /// rate ( /mm**2 /s)
    virtual double rate()const;
    
    /// set the static pointer 
    static void mgr(FluxMgr* );
    
    /// set the area of the target
    virtual void setTargetArea( double area);
    
    /// retrieve the area (a static, same for all fluxes)
    double targetArea()const;
    
    /// find which spectrum created the current particle
    virtual std::string findSource()const;
    
    /// return a unique number correcponding to that spectrum
    virtual int numSource()const;
    
    // virtual void addFactory( const IFactory* factory );
    
    virtual void addFactory(std::string name, const ISpectrumFactory* factory );
    /// get the transformtation matrix - the rest of these functions are now deprecated
    virtual CLHEP::HepRotation transformToGlast(double seconds, astro::GPS::CoordSystem index)const;
    
    
    //    insert(std::make_pair<std::string, const ISpectrumFactory*>(name, factory));
    
    EventSource* currentEvent(){return m_event;}
    EventSource* currentFlux(){return m_flux;}
    /// write the characteristics of the current source distribution to a stream
    void writeSourceCharacteristic(std::ostream& out);
    
    bool invalid()const;
private:
    
    EventSource* m_event;  
    EventSource* m_flux; // actual FluxSource used 

    static FluxMgr* s_mgr;
    
};

#endif
