/** 
 * @file LaunchDirection.h
 *
 * $Header$
 */

#ifndef _FluxSource_LaunchDirection_h
#define _FluxSource_LaunchDirection_h

#include "astro/SkyDir.h"

#include "astro/GPS.h"

#include <algorithm>
#include <sstream>

/** @class LaunchDirection
@brief nested launch strategy base class
$Header$
*/
class LaunchDirection  {
public:
    LaunchDirection():m_skydir(false),m_radius(0){}

    virtual ~LaunchDirection(){}

    LaunchDirection(double theta, double phi, std::string frame,double radius=0)
        :m_skydir(false)
        , m_radius(radius*M_PI/180),m_frame(frame)
    {
        HepVector3D dir(sin(theta)*cos(phi),sin(theta)*sin(phi),cos(theta));
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
    virtual void execute(double KE, double time){
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

    const HepVector3D& operator()()const {return dir();}

    virtual const HepVector3D& dir()const {
        static HepVector3D rdir;
        rdir = m_rottoglast * m_dir;
        if( m_radius>0 ) {
            // spread uniformly about a disk
            // rotate about perpendicular then about the original 
            HepVector3D t(rdir);
            t.rotate( m_radius*(sqrt(RandFlat::shoot())),  rdir.orthogonal()),  // rotate about the orthogonal
            t.rotate( RandFlat::shoot( 2*M_PI ), rdir); // rotate about the original direction
            rdir = t; //replace 
        }
        return rdir;
    }

    void setDir(const HepVector3D& dir){m_dir=dir;}


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
            astro::SkyDir zenDir(GPS::instance()->RAZenith(),GPS::instance()->DECZenith());
            return -m_dir*zenDir();
        }
        //if the direction is local
        return 1.0;
    }

    virtual const HepVector3D& skyDirection()const { return m_dir; }

private:
    HepRotation m_rottoglast;
    HepVector3D m_dir;
    bool  m_skydir;
    HepVector3D m_t;
    double m_radius;
    std::string m_frame;

};

#endif // _FluxSource_LaunchDirection_h
