/** 
 * @file LaunchPoint.h
 *
 * $Header$
 */

#ifndef _FluxSource_LaunchPoint_h
#define _FluxSource_LaunchPoint_h

#include "astro/SkyDir.h"

#include "astro/GPS.h"

/** @class LaunchPoint
@brief nested launch strategy base class for point determination

The virtual base class manages the point itself
$Header$
*/
class LaunchPoint  { 
public:
    LaunchPoint(){}
    LaunchPoint(const HepPoint3D& pt):m_pt(pt){}
    virtual ~LaunchPoint(){}

    /// access to direction, perhaps set by the execute()
    virtual const HepPoint3D& point()const {return m_pt;}
    const HepPoint3D& operator()()const{return point();}

    /// execute the strategy, perhaps depending on direction
    virtual void execute(const HepVector3D& ){};

    /// set the point
    void setPoint(const HepPoint3D& pt){ m_pt = pt;}

    /// return info, default if not overriden
    virtual std::string title()const{
        std::stringstream t;
        t << "point" << m_pt;
        return t.str();
    }

private:
    HepPoint3D m_pt;
};

#endif // _FluxSource_LaunchPoint_h
