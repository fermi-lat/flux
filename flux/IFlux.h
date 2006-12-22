// $Header: /nfs/slac/g/glast/ground/cvs/flux/flux/IFlux.h,v 1.7 2006/11/07 03:34:28 burnett Exp $

#ifndef _H_IFlux_
#define _H_IFlux_

// includes
#include <string>
#include "CLHEP/Vector/Rotation.h"
#include "CLHEP/Vector/ThreeVector.h"
#include "astro/GPS.h"


class ParticleProperty;
class EventSource;

class ISpectrumFactory;

/** 
* \class IFlux
* \brief The virtual interface for Flux-type objects.
*
* \author Toby Burnett tburnett@u.washington.edu
* 
  Abstract interface for an object that generates particles, Flux

  * $Header: /nfs/slac/g/glast/ground/cvs/flux/flux/IFlux.h,v 1.7 2006/11/07 03:34:28 burnett Exp $
*/
class IFlux {
public:
    /// ctor, select the name
    IFlux(std::string =""){};
    virtual ~IFlux(){}
    
    /// name of the flux
    virtual std::string name()const=0;
    
    /// full title of the flux
    virtual std::string title()const = 0;
    
    /// generate a new entry trajectory. Must be valid, that is, there are available sources. 
    /// May become invalid as a result of this call, return is validity
    /// if invalid, std::runtime exception will be thrown
    virtual bool generate()=0;

    /// the particle name of the last particle generated 
    virtual std::string particleName()const=0;
    
    /// the particle property entry for the last particle generated 
    //virtual ParticleProperty* property()const=0;
    
    /// its kinetic energy
    virtual double energy()const=0;
    
    /// starting point 
    virtual Hep3Vector launchPoint()const=0;
    
    /// direction
    virtual Hep3Vector launchDir()const=0;
    
    /// time (s) (absolute or elapsed??)
    virtual double time()const=0;
    
    /// return rate ( /mm**2 /s)
    virtual double rate()const=0;
    
    /// set the area of the target
    virtual void setTargetArea( double area)=0;
    
    /// get the target area
    virtual double targetArea()const =0;
    
    /// find which spectrum created the current particle
    virtual std::string findSource()const=0;
    
    /// return a unique number correcponding to that spectrum
    virtual int numSource()const=0;
    
    /// pass a specific amount of time
    virtual void pass ( double t)=0;
        
    virtual void addFactory(std::string name, const ISpectrumFactory* factory )=0;
    
    virtual /*int*/double gpsTime()const=0;
    
    virtual EventSource* currentEvent()=0;
    virtual EventSource* currentFlux()=0;

    /// write the characteristics of the current source distribution to a stream
    virtual void writeSourceCharacteristic(std::ostream& out)=0;

    /// get the transformtation matrix - the rest of these functions are now deprecated
    virtual CLHEP::HepRotation transformToGlast(double seconds,astro::GPS::CoordSystem index)const=0;
    
    /// may become invalid if no sources are enabled
    virtual bool invalid()const=0; 

};


#endif // _H_FluxSvc
