/** @file CompositeSource.cxx
@brief Define CompositeSource

$Header: /nfs/slac/g/glast/ground/cvs/flux/src/CompositeSource.cxx,v 1.18 2008/01/07 22:58:06 burnett Exp $
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
, m_recent(0),m_occulted(false)
{
}

CompositeSource::~CompositeSource()
{
    for (std::vector<EventSource*>::iterator it = m_sourceList.begin();
        it != m_sourceList.end(); ++it ) delete (*it);
}


void CompositeSource::addSource (EventSource* aSource)
{
    int index( m_sourceList.size() ); 
    m_sourceList.push_back(aSource);
    EventSource* dummy(0);
    // insert in the map, tagged as needing to be evaluated
    m_source_map.insert(std::make_pair( -1, std::make_pair(aSource, dummy)));

    // default identifier is the index
    m_ident[aSource] = index;

}

EventSource* CompositeSource::process(CompositeSource::SourceMap::iterator it, double time)
{
    // set up defaults to flag that this source was used, and will be rerun on next iteration
    EventSource * source(0);     // will set to the actual source if updating
    double nexttime(-1);         // will be first 

    // now get details and remove this entry from the map
    m_recent = it->second.first;  // the member
    EventSource* actual = it->second.second;   // the actual FluxSource object if member is Composite
    m_source_map.erase(it);       // remove the entry from the map

    if( actual==0){
        // node needs to be run, so ask for next photon
        source = m_recent->event(time);
        if( ! m_recent->enabled() ) return 0; // if disabled, no more calls

        double nextinterval( m_recent->interval() );
        nexttime = time+nextinterval; 
        source->setTime(nexttime); // tell event its time
    }
    m_source_map.insert(std::make_pair( nexttime, std::make_pair(m_recent, source)));
    return actual; // now return either a valid source, or null

}

EventSource* CompositeSource::event (double time)
{

    if( !enabled()){
        throw std::runtime_error("CompositeSource::event called when disabled");
    }
    
    EventSource* actual(0);
    do {
        // get earliest entry (
        SourceMap::iterator it= m_source_map.begin();
        if( it==m_source_map.end()){ // no sources left
            disable();
            return this;
        }
        // process it, returning either null or a valid source
        actual = process(it, time); 

    }while (actual==0); // loop until valid source

    double nexttime( actual->time());
    if( nexttime < time){
        throw std::runtime_error("CompositeSource::event: invalid time");
    }

    // save the actual interval and time and return the current actual source
    setInterval(nexttime-time);
    setTime(nexttime);
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
    // get the value from the m_ident map, using the most recent source as a key
    std::map<EventSource*, unsigned int>::const_iterator it = m_ident.find(m_recent);
    int index( it->second );
    // now call this function associated with that source:
    int t( m_recent->numSource() ); 
    // this will be -1 if the source is not Composite. if it is composite, it will be
    // 1000 times its index
    if( t!= -1){
        // it is composite. Get its recent guy, 
        EventSource* actual = dynamic_cast<CompositeSource*>(m_recent)->m_recent;
        if( actual !=0) { // check: maybe 
            int id (actual->identifier());
            if( id>0 ) {
                return id;
            }
        }
    }
    return EventSource::s_id_offset + 1000*index + (t==-1? 0:  t/1000);

}

