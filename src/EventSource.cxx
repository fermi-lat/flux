/** @file EventSource.cxx
    @brief Implementation of class EventSource

   $Header: /nfs/slac/g/glast/ground/cvs/FluxSvc/src/EventSource.cxx,v 1.17 2003/03/21 19:14:37 jrb Exp $
*/

#include "flux/EventSource.h"

#include <xercesc/dom/DOM_Element.hpp>
#include "xml/Dom.h"
#include "flux/FluxException.h"

#include <sstream>

unsigned int  EventSource::s_id = 0;
double  EventSource::s_total_area = 6.; // area in m^2

EventSource::EventSource (double aFlux, unsigned acode)
:  m_enabled(true), m_flux(aFlux),  m_code(acode)
{
    std::stringstream  s;
    
    s << "Source_" << (++s_id) << '\0';
    if (acode == 0) code(s_id); // automatically assign event codes...
    
    m_name = s.str();
}
EventSource::~EventSource()
{}

double EventSource::flux (double time) const
{
  // Purpose and Method: This method returns the flux of the particular source.
  // Inputs  - current time
  // Outputs - flux, in units of (particles/(m^2*sr*sec))
    return m_flux;  // default if not overridden
}


double  EventSource::rate (double time )const
{
  // Purpose and Method: This method returns the rate of particles entering the detector.
  // Inputs  - current time
  // Outputs - rate, in units of (particles/sec)
    return enabled()? (solidAngle()*flux(time)*s_total_area) :0;
}


double	EventSource::solidAngle () const{
    return m_solid_angle;
}

// UI titles
std::string EventSource::fullTitle () const 
{ return std::string("EventSource");   }
std::string EventSource::displayTitle () const  {  return m_name; }

// inline function declarations:


const std::string& EventSource::name () const	{   return m_name;  }
void EventSource::name (const std::string& value)    { m_name = value;   }

double    EventSource::totalArea () { return s_total_area; }
void    EventSource::totalArea (double value) { s_total_area = value; }

unsigned EventSource::code () const { return m_code; }
void EventSource::code ( unsigned c ) { m_code = c; }

void EventSource::setFlux(double value){ m_flux=value; }
