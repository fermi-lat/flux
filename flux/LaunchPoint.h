/** 
 * @file LaunchPoint.h
 * @brief Declare LaunchPoint class
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/flux/flux/LaunchPoint.h,v 1.3 2006/03/21 01:28:55 usher Exp $
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
    LaunchPoint(const CLHEP::Hep3Vector& pt):m_pt(pt){}
    virtual ~LaunchPoint(){}

    /// access to direction, perhaps set by the execute()
    virtual const CLHEP::Hep3Vector& point()const {return m_pt;}
    const CLHEP::Hep3Vector& operator()()const{return point();}

    /// execute the strategy, perhaps depending on direction
    virtual void execute(const CLHEP::Hep3Vector& ){};

    /// set the point
    void setPoint(const CLHEP::Hep3Vector& pt){ m_pt = pt;}

    /// return info, default if not overriden
    virtual std::string title()const{
        std::stringstream t;
        t << "point" << m_pt;
        return t.str();
    }

private:
    CLHEP::Hep3Vector m_pt; ///< the point we manage
};

#endif // _FluxSource_LaunchPoint_h
