/** 
 * @file LaunchDirection.h
 * @brief Declare LaunchDirection class
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/flux/flux/LaunchDirection.h,v 1.5 2006/03/21 01:28:55 usher Exp $
 */

#ifndef _FluxSource_LaunchDirection_h
#define _FluxSource_LaunchDirection_h

#include "CLHEP/Geometry/Point3D.h"
#include "CLHEP/Geometry/Vector3D.h"
#include "CLHEP/Random/RandFlat.h"
#include "CLHEP/Vector/Rotation.h"

// Hack for CLHEP 1.9.2.2
#ifndef HepVector3D
namespace HepGeom {
    typedef Vector3D<double> HepVector3D;
    typedef Point3D<double>  HepPoint3D;
}
#endif

#include "astro/SkyDir.h"

#include "astro/GPS.h"

#include <algorithm>
#include <sstream>

/** @class LaunchDirection
@brief launch strategy base class

*/
class LaunchDirection  {
public:
    LaunchDirection():m_skydir(false),m_radius(0){}

    virtual ~LaunchDirection(){}

    LaunchDirection(double theta, double phi, std::string frame,double radius=0)
        :m_skydir(false)
        , m_radius(radius*M_PI/180),m_frame(frame)
    {
        HepGeom::HepVector3D dir(sin(theta)*cos(phi),sin(theta)*sin(phi),cos(theta));
        setDir(-dir); // minus due to z axis pointing UP!
    }
    LaunchDirection(astro::SkyDir sky, double radius=0)
        :m_skydir(true)
        , m_radius(radius*M_PI/180)
    {
        m_dir = -sky.dir();
    }
    /** @brief choose a direction
    @param KE kinetic energy
    @param time mission time
    */
   virtual void execute(double /*KE*/, double time){
       using astro::GPS;
        if(m_skydir){
            //here, we have a SkyDir, so we need the transformation from a SkyDir to GLAST.
            m_rottoglast = GPS::instance()->transformToGlast(time,GPS::CELESTIAL);//->transformCelToGlast(time);
        }else{
            if(m_frame=="zenith"){
                //The direction is in the earth zenith system, and the rotation to GLAST is needed:
                m_rottoglast = GPS::instance()->transformToGlast(time,GPS::ZENITH);
            }else{
                //otherwise, the direction is in the spacecraft system, and the rotation to GLAST is the identity:
                m_rottoglast = GPS::instance()->transformToGlast(time,GPS::GLAST);
            }
        }
    }

    const HepGeom::HepVector3D& operator()()const {return dir();}

    virtual const HepGeom::HepVector3D& dir()const {
        static HepGeom::HepVector3D rdir;
        rdir = m_rottoglast * m_dir;
        if( m_radius>0 ) {
            // spread uniformly about a disk
            // rotate about perpendicular then about the original 
            HepGeom::HepVector3D t(rdir);
            t.rotate( m_radius*(sqrt(CLHEP::RandFlat::shoot())),  rdir.orthogonal()),  // rotate about the orthogonal
            t.rotate( CLHEP::RandFlat::shoot( 2*M_PI ), rdir); // rotate about the original direction
            rdir = t; //replace 
        }
        return rdir;
    }

    void setDir(const HepGeom::HepVector3D& dir){
        m_dir=dir;
    }


    //! solid angle: default of 1. for a point source
    virtual double solidAngle()const {
        return 1.;
    }

    /// return info, default if not overriden
    virtual std::string title()const{
        std::stringstream t;
        t << " dir" << m_dir ;
        if( m_radius>0 ) { t << " radius " << m_radius ;}
        return t.str();
    }

    /// return the cosine of the angle between the incoming direction and the earth's zenith
    virtual double zenithCosine()const{
        if(m_skydir){
            astro::SkyDir zenDir(astro::GPS::instance()->RAZenith(),astro::GPS::instance()->DECZenith());
            return -m_dir*zenDir();
        }
        //if the direction is local
        return 1.0;
    }

    virtual const HepGeom::HepVector3D& skyDirection()const { return m_dir; }

private:
    CLHEP::HepRotation m_rottoglast;
    HepGeom::HepVector3D m_dir;
    bool  m_skydir;
    HepGeom::HepVector3D m_t;
    double m_radius;
    std::string m_frame;

};

#endif // _FluxSource_LaunchDirection_h
