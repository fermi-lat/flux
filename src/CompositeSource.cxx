/** @file CompositeSource.cxx
@brief Define CompositeSource

$Header: /nfs/slac/g/glast/ground/cvs/flux/src/CompositeSource.cxx,v 1.14 2008/01/06 22:04:30 burnett Exp $
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

void CompositeSource::map_insert(double time, EventSource* member, EventSource * source=0)
{
    m_source_map.insert(std::make_pair( time, std::make_pair(member, source)));
}
void CompositeSource::addSource (EventSource* aSource)
{
    m_sourceList.push_back(aSource);
#if 1 // new map-based code
    // insert in the map, tagged as needing to be evaluated
    map_insert(-1.,  aSource, 0);
    // tag identifier
    m_ident[aSource] = m_sourceList.size()-1;

#else
    //here, set up the associated vectors by default.
    m_unusedSource.push_back(0);
    m_sourceTime.push_back(-1);
    m_eventList.push_back(0);
#endif
}

EventSource* CompositeSource::event (double time)
{

    if( !enabled()){
        throw std::runtime_error("CompositeSource::event called when disabled");
    }
#if 1 // new map-based code
    EventSource* actual(0);
    m_recent=0;
    double nexttime(0);
    do {
        // get next one, and remove from the map
        SourceMap::iterator it= m_source_map.begin();
        if( it==m_source_map.end()){
            // no sources left
            disable();
            return this;
        }
        nexttime = it->first; if(nexttime<0) nexttime=time; // initial
        m_recent = it->second.first;
        actual = it->second.second;
        m_source_map.erase(it);

        // evaluate its next event time, put back into queue
        EventSource * source = m_recent->event(nexttime);
        if( m_recent->enabled()){
            double newtime(nexttime+m_recent->interval() ); 
            map_insert(newtime, m_recent, source);
        }

    }while (actual ==0 );// loop until first evaluation of all sources

    if( nexttime < time){
        throw std::runtime_error("CompositeSource::event: invalid time");
    }
    // save the actual interval and return the current source
    setInterval(nexttime-time);
    m_numofiters = m_ident[m_recent]; // the sequence of this guy

    m_occulted=actual->occulted();

    return actual;
#else // old code
    int i=0; //for iterating through the m_unusedSource vector
    int winningsourcenum=-1; //the number of the "winning" source
    EventSource::setTime(time);

    m_numofiters=0;

    // more than one:: choose on basis of relative rates
    std::vector<EventSource*>::iterator  now = m_sourceList.begin();
    std::vector<EventSource*>::iterator  it = now;

    double intrval=0.,intrmin=1e10; //100000.;
    int q;
    for (q=0 ; now != m_sourceList.end(); ++now, ++i, ++q) {
        if( ! (*now)->enabled() ) continue; // ignore if turned off
        if(m_unusedSource[i]==1){
            // was not used yet: update the interval
            intrval=m_sourceTime[i]-time;
        }else{

            // this was was used last time: get a new interval, unless now disabled

            EventSource* candidate = (*now)->event(time);
            if( !candidate->enabled() ) continue; // skip this guy, no longer active
            m_eventList[i] = candidate; // to initialize particles, so that the real interval for the particle is gotten.
            intrval=m_sourceList[i]->interval(time);
            if( intrval <=0 ){
                throw("CompositeSource::event: zero or negative interval");
            }
            m_unusedSource[i]=1;
            m_sourceTime[i]=time + intrval;
        }

        if(intrval < intrmin){
            //the present source is "winning" here
            it=now;
            intrmin=intrval;
            m_numofiters=q;
            winningsourcenum=i;
        }

        m_recent = (*it);
    }
    if( winningsourcenum<0 ) {
        // nothing left: we are disabled.
        disable();
        return this;
    }
    //note:the internal interval() function takes absolute time.
    setInterval(time+intrmin);
    m_unusedSource[winningsourcenum]=0; //the current "winning" source is getting used...
    //now, check to see if this source is occulted:
    m_occulted=m_eventList[winningsourcenum]->occulted();
    // now ask the chosen one to return the event.
    return m_eventList[winningsourcenum];
#endif
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

