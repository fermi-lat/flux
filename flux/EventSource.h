/** @file EventSource.h
   @brief Declaration of EventSource

   $Header: /nfs/slac/g/glast/ground/cvs/flux/flux/EventSource.h,v 1.2 2003/10/01 22:21:50 srobinsn Exp $
*/

#ifndef EventSource_h
#define EventSource_h 1
#include "CLHEP/Geometry/Point3D.h"
#include "CLHEP/Geometry/Vector3D.h"


/** 
* \class EventSource
*
* \brief  Base class for managing  sources.

This the abstract base class for source, (FluxSource) or a list of sources (CompositeSource)
* 
* $Header: /nfs/slac/g/glast/ground/cvs/flux/flux/EventSource.h,v 1.2 2003/10/01 22:21:50 srobinsn Exp $
*/


#include <string>

class EventSource
{
public:
    /// ctor/dtor
    EventSource (double aFlux = 1.0, unsigned acode = 0);
    virtual ~EventSource();
    
    ///    a randomized interval to the next event - default is 1/rate()
    virtual double interval (double) = 0;
    
    ///    calculate the rate for a given flux/solid angle integral (NOTE: integral of solid angle)
    // virtual double  rate ( double solid_angle, double flux );	
    virtual double  rate (double time)const;

    ///    abstract method - create an event
    virtual EventSource* event (double) = 0;	

    ///    UI titles - used for tuple header (verbose) or window title (display)
    virtual std::string fullTitle () const;
    virtual std::string displayTitle () const;
    
    ///    flux for this source in (p/(m^2*sr*sec))
    virtual double	flux (double time) const;
    virtual void      setFlux (double value);
    
    ///    disable/enable, test this particular source 
    void      disable (){m_enabled=false;}
    void      enable(){m_enabled=true;}
    bool enabled()const{return m_enabled;}
    
    ///    integral of solid angle over which flux is incident
    virtual double	solidAngle () const;
    
    ///	name of this flux source - for UI
    const std::string& name () const;
    void name (const std::string& value);
    
    ///    code - for monte-carlo study
    unsigned  code () const;
    virtual void  code ( unsigned );
    
    ///    area 
    static double	totalArea ();
    static void	totalArea ( double value );
	/// Is the current particle occulted?
	virtual bool occulted() = 0;
    /// virtual event number: should be filled in by subclass
    virtual int eventNumber()const{return -1;} 
    
    ///  say which source created the current particle
    virtual std::string findSource()const{return "";}
    
    
    /// return a unique number correcponding to that spectrum
    virtual int  numSource()const{return -1;}
    
    virtual double time()const{return m_time;}
    virtual void setTime(double time){m_time=time;}
    
    
    //return how many sources are in the sourcelist (defaults to 1 if only a single FluxSource)
    virtual int howManySources(){return 1;}

    /// write the characteristics of the current source distribution to a stream
    virtual void writeSourceCharacteristic(std::ostream& out){
        out << fullTitle() << std::endl;
        out<< "default message - no sources" << std::endl;
    }
    /// the current particle:: expect to be overridden by FluxSource
    virtual std::string particleName() { return std::string("unknown");}

    virtual double energy()const { return 0;}
    virtual const HepVector3D& launchDir()const { static HepVector3D dummy; return dummy;}
    virtual const HepPoint3D&  launchPoint()const { static HepPoint3D dummy; return dummy;}

    
private:
    double m_time;    // elapsed time, really only needed for EventSource
    
    bool m_enabled;           // toggle that it is enabled
    double m_flux;		// representative flux for this event source...
    std::string m_name;       // name of the event source (UI)
    unsigned  m_code;         // code id for this event source - for identification in the tuple.
    
    static unsigned int  s_id;    // id for new EventSources...
    static double s_total_area;   // total area for flux generation (in square meters)
    double m_solid_angle;
};


#endif
