/** @file EventSource.cxx
    @brief Implementation of class EventSource

   $Header: /nfs/slac/g/glast/ground/cvs/ScienceTools-scons/flux/src/EventSource.cxx,v 1.10 2008/01/14 20:16:07 burnett Exp $
*/

#include "flux/EventSource.h"

#include "flux/FluxException.h"

#include <sstream>

int  EventSource::s_id = 0;
double  EventSource::s_total_area = 6.; // area in m^2
double  EventSource::s_backoff = 2000.; // in mm
int EventSource::s_id_offset=0;     // offset for assigning ids in CompositeSource

bool EventSource::s_applyAlign(false);
CLHEP::HepRotation EventSource::s_alignMatrix;



std::vector<double> EventSource::s_cone; 


EventSource::EventSource (double aFlux, unsigned acode)
   :  m_enabled(true), m_flux(aFlux),  m_code(acode), m_applyEdisp(true)
{
    std::stringstream  s;
    
    s << "Source_" << (++s_id) << '\0';
    if (acode == 0) code(s_id); // automatically assign event codes...
    
    m_name = s.str();
}
EventSource::~EventSource()
{}

auto EventSource::flux (double time) const -> double
{
  // Purpose and Method: This method returns the flux of the particular source.
  // Inputs  - current time
  // Outputs - flux, in units of (particles/(m^2*sr*sec))
    return m_flux;  // default if not overridden
}

auto EventSource::interval()const -> double
{
    return m_interval;
}

auto EventSource::setInterval (double time) -> double
{
    return (m_interval = time);
    if( time<=0 ){
        std::cout << "interval set <=0: " << time << std::endl;
    }
}


auto  EventSource::rate (double time )const -> double
{
  // Purpose and Method: This method returns the rate of particles entering the detector.
  // Inputs  - current time
  // Outputs - rate, in units of (particles/sec)
    return enabled()? (solidAngle()*flux(time)*s_total_area) :0;
}


auto	EventSource::solidAngle () const -> double{
    return m_solid_angle;
}

// UI titles
auto EventSource::fullTitle () const -> std::string 
{ return std::string("EventSource");   }
auto EventSource::displayTitle () const -> std::string  {  return m_name; }

// inline function declarations:


auto EventSource::name () const -> std::string	{   return m_name;  }

void EventSource::setName (const std::string& value)    { m_name = value;   }

auto    EventSource::totalArea () -> double { return s_total_area; }
void    EventSource::totalArea (double value) { s_total_area = value; }

auto EventSource::code () const -> int { 
    return m_code; 
}
void EventSource::code ( int c ) { 
    m_code = c; 
}

void EventSource::setFlux(double value){ m_flux=value; }

void EventSource::setAlignmentRotation(const CLHEP::HepRotation &align)
{
    s_alignMatrix = align;
    s_applyAlign = true;
}


