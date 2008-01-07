/** @file CompositeSource.cxx
@brief Define CompositeSource

$Header: /nfs/slac/g/glast/ground/cvs/flux/src/CompositeSource.cxx,v 1.16 2008/01/07 12:14:28 burnett Exp $
*/

#include "flux/CompositeSource.h"  


#include <sstream>
#include <cassert>
#include <numeric> // for accumulate
#include <functional>
#include <iomanip>
#include <cmath>
#include <stdexcept>
#include <cassert>

CompositeSource::CompositeSource (double aRate)
: EventSource(aRate)
, m_numofiters(0)
, m_recent(0),m_occulted(false)
{
}

CompositeSource::~CompositeSource()
{
    for (std::vector<EventSource*>::iterator it = m_sourceList.begin();
        it != m_sourceList.end(); ++it ) delete (*it);
}

void CompositeSource::map_insert(double time, EventSource* member, EventSource*actual)
{
    EventSource * source(0); // will set to the actual source if updating
    double nexttime(time);
    if( actual==0){
        source = member->event(time);
        if( member->enabled()){
            double nextinterval( member->interval() );
            nexttime = time+nextinterval; 
            source->setTime(nexttime); // tell event its time
        }
    }else{
        nexttime=-1;
    }
    m_source_map.insert(std::make_pair( nexttime, std::make_pair(member, source)));
}
void CompositeSource::addSource (EventSource* aSource)
{
    m_sourceList.push_back(aSource);
    // insert in the map, tagged as needing to be evaluated
    map_insert(-1.,  aSource);
    // tag identifier
    m_ident[aSource] = m_sourceList.size()-1;

}

EventSource* CompositeSource::event (double time)
{

    if( !enabled()){
        throw std::runtime_error("CompositeSource::event called when disabled");
    }
    
    EventSource* actual(0);
    double nexttime(0);
    do {
        // get most recent entry, and remove from the map
        SourceMap::iterator it= m_source_map.begin();
        if( it==m_source_map.end()){
            // no sources left
            disable();
            return this;
        }
        nexttime = it->first;         // initial
        m_recent = it->second.first;  // the member
        actual = it->second.second;   // the actual FluxSource object if member is Composite
        m_source_map.erase(it);
        map_insert( time, m_recent, actual);    // will catch on next iteration
        
    }while (nexttime< 0); // loop until evaluation of all sources

    if( nexttime < time){
        throw std::runtime_error("CompositeSource::event: invalid time");
    }
    // save the actual interval and return the current source
    setInterval(nexttime-time);
    m_numofiters = m_ident[m_recent]; // the sequence of this guy

    m_occulted=actual->occulted();

    return actual;

}

std::string CompositeSource::fullTitle () const
{
    std::stringstream  s;
    std::vector<EventSource*>::const_iterator	it = m_sourceList.begin();

    while (it != m_sourceList.end()) {

        s << (*it)->fullTitle() << " ";
        ++it;
        if (it != m_sourceList.end())    s << "+ ";
    }
    std::string t(s.str());
    return t;
}

std::string CompositeSource::displayTitle () const
{
    return (m_recent == 0) ? "" : m_recent->displayTitle();
}

double CompositeSource::rate(double time) const
{
    //m_time += m_time-time;
    std::vector<EventSource*>::const_iterator it = m_sourceList.begin();
    double	total_rate = 0.;
    for(;it != m_sourceList.end();++it) {
        double rr = fabs((*it)->rate(time));
        total_rate += rr;
    }
    return total_rate;
}

void CompositeSource::printOn(std::ostream& out)const
{
    out << "Source(s), total rate="<< rate(EventSource::time()) << std::endl;

    for( std::vector<EventSource*>::const_iterator it = m_sourceList.begin();
        it != m_sourceList.end();++it)	{
            out <<  std::setw(8) << std::setprecision(4) << (*it)->rate(EventSource::time()) <<" Hz, "
                << '#' << std::setw(6) << (*it)->eventNumber() <<' '
                << (*it)->name() << ' '<< (*it)->fullTitle() << std::endl;

        }
}

std::string CompositeSource::findSource()const
{
    return m_recent->fullTitle();
}

int  CompositeSource::numSource()const
{
    ///Purpose: Return a unique number correcponding to the current spectrum.
    // if selected source is composite itself, (id=-1)  add its id as an integer
    int t=m_recent->numSource();
    return EventSource::s_id_offset + 1000*m_numofiters + (t==-1? 0:  t/1000);
}

