/** @file CompositeSource.h
    @brief CompositeSource declaration
    
  $Header: /nfs/slac/g/glast/ground/cvs/flux/flux/CompositeSource.h,v 1.2 2003/10/01 22:21:50 srobinsn Exp $
*/

#ifndef CompositeSource_h
#define CompositeSource_h 1
/** 
* \class CompositeSource
*
* \brief  holds multiple Eventsource objects ; acts as a container for them.
* Each time an event() is called, CompositeSource goes through a process of deciding 
* "which source" it is representing this time.  Old particles are held, along with the
* time of their arrival, until use.
* 
* $Header: /nfs/slac/g/glast/ground/cvs/flux/flux/CompositeSource.h,v 1.2 2003/10/01 22:21:50 srobinsn Exp $
*/

#include "flux/EventSource.h"
#include <vector>

class CompositeSource : public EventSource { 
public:
    ///    constructor/destructor
    CompositeSource (double aRate = 1.0);
    virtual ~CompositeSource();
    
    
    ///    add a source to the list
    virtual void addSource (EventSource* aSource);
    
    /// generate an event from from one of the sources 
    /// which make up the composite, and return a pointer to it
    virtual EventSource* event (double time);
    
    /// rate - compute overall rate...
    virtual double rate (double time)const;
    
    /// flux into 1 m^2 integrated over angles
    virtual double flux(double time)const{
        return rate(time)/totalArea();}
    
    ///    full-length title description of this EventSource.
    virtual std::string fullTitle () const;
    
    ///    brief title description (for display) for this event source.
    virtual std::string displayTitle () const;
    
    /// dump current list of sources, rates
    void printOn(std::ostream& out)const;
    
    ///  say which source created the current particle
    std::string findSource()const;
    
    /// return a unique number correcponding to that spectrum
    virtual  int numSource()const;
    
    
    ///	    list of sources which make up this composite
    std::vector< EventSource* >& sourceList ();
    const std::vector< EventSource* >& sourceList () const;
    void sourceList (const std::vector< EventSource* >& value);
    
    /// interval to the next event
    double interval (double time){return m_interval - time;}
    
    /// set the absolute time to the next event
    double setInterval (double time){return (m_interval = time);}
    
    /// double m_time; 
    double m_interval;
    
    /// return how many sources are in the sourcelist
    int howManySources(){return m_sourceList.size();}

    /// is the most recent photon occulted?
    bool occulted(){return m_occulted;}
    
protected:
    
    //number of times we've iterated the front() pointer into sourcelist 
    //to get the current particle - represents the source
    int m_numofiters;
    
    //private: 
    std::vector< EventSource* > m_sourceList;
    std::vector< EventSource* > m_eventList;

    //vector of flags, holds whether or not the current source has a remaining unused particle.
    std::vector<int> m_unusedSource;
    //vector of recorded arrival times of held sources.
    std::vector<double> m_sourceTime;
    EventSource*  m_recent;
	//is the photon from the most recent source occulted?
	bool m_occulted;
};

inline std::vector< EventSource* >& CompositeSource::sourceList ()
{
    return m_sourceList;
}

inline const std::vector< EventSource* >& CompositeSource::sourceList () const
{
    return m_sourceList;
}

inline void CompositeSource::sourceList (const std::vector< EventSource* >& value)
{
    m_sourceList = value;
}
#endif
