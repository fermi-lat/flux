/** 
 * @file LaunchPoint.h
 * @brief Declare LaunchPoint class
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/flux/flux/LaunchPoint.h,v 1.2 2005/05/05 16:50:06 burnett Exp $
 */

#ifndef _FluxSource_LaunchPoint_h
#define _FluxSource_LaunchPoint_h

#include "astro/SkyDir.h"

#include "astro/GPS.h"

/** @class LaunchPoint
@brief launch strategy base class for point determination

The virtual base class manages the point itself
*/
class LaunchPoint  { 
public:
    LaunchPoint(){}
    LaunchPoint(const HepGeom::HepPoint3D& pt):m_pt(pt){}
    virtual ~LaunchPoint(){}

    /// access to direction, perhaps set by the execute()
    virtual const HepGeom::HepPoint3D& point()const {return m_pt;}
    const HepGeom::HepPoint3D& operator()()const{return point();}

    /// execute the strategy, perhaps depending on direction
    virtual void execute(const HepGeom::HepVector3D& ){};

    /// set the point
    void setPoint(const HepGeom::HepPoint3D& pt){ m_pt = pt;}

    /// return info, default if not overriden
    virtual std::string title()const{
        std::stringstream t;
        t << "point" << m_pt;
        return t.str();
    }

private:
    HepGeom::HepPoint3D m_pt; ///< the point we manage
};

#endif // _FluxSource_LaunchPoint_h
