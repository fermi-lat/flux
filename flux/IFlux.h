// $Header: /nfs/slac/g/glast/ground/cvs/FluxSvc/FluxSvc/IFlux.h,v 1.28 2003/03/20 19:55:32 burnett Exp $

#ifndef _H_IFlux_
#define _H_IFlux_

// includes
#include <string>
#include "CLHEP/Geometry/Point3D.h"
#include "CLHEP/Geometry/Vector3D.h"
#include "CLHEP/Vector/Rotation.h"

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

  * $Header: /nfs/slac/g/glast/ground/cvs/FluxSvc/FluxSvc/IFlux.h,v 1.28 2003/03/20 19:55:32 burnett Exp $
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
    
    /// generate a new entry trajectory
    virtual void generate()=0;
    
    /// the particle name of the last particle generated 
    virtual std::string particleName()const=0;
    
    /// the particle property entry for the last particle generated 
    //virtual ParticleProperty* property()const=0;
    
    /// its kinetic energy
    virtual double energy()const=0;
    
    /// starting point 
    virtual HepPoint3D launchPoint()const=0;
    
    /// direction
    virtual HepVector3D launchDir()const=0;
    
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
    
    ///get the transformation matrix due to orientation of the Galaxy 
    virtual HepRotation CELTransform(double time)const=0;
    
    ///get the transformation matrix due to orientation of the spacecraft.
    virtual HepRotation orientTransform(double time)const=0;
    
    virtual void addFactory(std::string name, const ISpectrumFactory* factory )=0;
    
    virtual /*int*/double gpsTime()const=0;
    
    ///this transforms glast-local (cartesian) vectors into galactic (cartesian) vectors
    virtual HepRotation transformGlastToGalactic(double time)const=0;
    
    virtual EventSource* currentEvent()=0;
    virtual EventSource* currentFlux()=0;

    /// write the characteristics of the current source distribution to a stream
    virtual void writeSourceCharacteristic(std::ostream& out)=0;
    
};


#endif // _H_FluxSvc
