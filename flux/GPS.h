// $Header: /nfs/slac/g/glast/ground/cvs/flux/flux/GPS.h,v 1.1.1.1 2003/07/29 18:22:14 burnett Exp $

#if !defined(_H_GPS_CLASS)
#define _H_GPS_CLASS



#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "facilities/Scheduler.h"
#include "facilities/Observer.h"

#include "astro/SkyDir.h"
#include "astro/EarthOrbit.h"
#include "CLHEP/Vector/Rotation.h"

#include <iostream>


/** 
* \class GPS
* \brief Models the Global Positoning System for a spacecraft. Handles time, position, and orientation for the instrument as a whole.
* 
* $Header: /nfs/slac/g/glast/ground/cvs/flux/flux/GPS.h,v 1.1.1.1 2003/07/29 18:22:14 burnett Exp $
 Represents the Global Positioning System on-board the spacecraft. An Orbit
  object is used to compute the spacecraft's position and pointing characteristics.
Time is tracked through this object, and synchronized with the Scheduler for 
discrete event simulation. An expansion factor is provided to allow for acceleration
or randomization of the orbit. If the expansion factor is negative, then the position
of the spacecraft is chosen as a random distribution over an orbit. Otherwise, the expansion
factor represents an acceleration of the spacecraft's orbit. Ie. an expansion factor
of 2 would reduce the orbit period of the spacecraft by 1/2.

*/
class GPS  
{
public:

    enum RockType { 
            NONE,  //!  No rocking rotation done at all.
            UPDOWN, //! Satellite will be rocked toward the north pole in the northern hemisphere, opposite in the south.
            SLEWING, //! (experimental) like UPDOWN, except that rotation at equator happens gradually.
            ONEPERORBIT, //! LAT rocked northward for one orbit, southward for the next.
            EXPLICIT, //!  Explicit angles given - this is used only if rotAngles get set.
			POINT, //!  The Lat points in a given direction.  Here, m_rotangles is in (l,b) format for pointing.
			HISTORY //! This setting is for using a previously generated pointing database to represent the orbit.
        };

    class Coords {
    public:
        Coords( double alat, double alon, double apitch
            , double ayaw, double aroll, GPStime atime, double phase ); 
        Coords();
        
        GPStime time () const { return m_time; }
        double  lat () const { return m_lat; }
        double  lon () const { return m_lon; }
        double  phase () const { return m_phase; }
        
    private:
        GPStime m_time;
        double m_lat, m_lon, m_pitch, m_yaw, m_roll, m_phase;
    };

    // const access
    
    /// GPS synchronized time for the satellite
    GPStime	time () const; 
    /// present latitude
    double	lat () const; 
    /// present longitude
    double	lon () const;  
    /// expansion of the current orbit
    double      expansion () const; 
    /// sample interval for random orbit distribution
    GPStime     sampleintvl () const; 
    /// return the orbit's ascending longitude 
    double     ascendingLon()const; 
    /// access m_rotangles
    std::pair<double,double> rotateAngles(); 

 
    // set data
    
    /// get the pointing characteristics of the satellite, given a location and rocking angle.
    void getPointingCharacteristics(double seconds);
    
    /// pass a specific amount of time
    void    pass ( double );
    /// set the expansion factor for the orbit (-1) = random
    void    expansion ( double );
    /// synchronize w. scheduler
    void    synch ();    
    /// set the sample interval
    void    sampleintvl ( GPStime );
    /// special to set the ascending longitude	
    void    ascendingLon(double);   
    /// set m_rotangles
    void    rotateAngles(std::pair<double,double> coords); 
    
    /// print time & position
    void    printOn(std::ostream& out) const; 
    
    // notification
    void notifyObservers() { m_notification.notify();}
    Subject&    notification() { return m_notification; }
    
    // static access/destruction
    static GPS*	instance();
    static void     kill ();
    
    /// return the rotation for compensation for the rocking angles.
    HepRotation rockingAngleTransform(double seconds);
    
    ///this transforms glast-local (cartesian) vectors into galactic (cartesian) vectors
    HepRotation transformGlastToGalactic(double seconds);

    HepRotation GPS::CELTransform(double seconds);

    HepRotation GPS::transformCelToGlast(double seconds);

    double rockingDegrees(double rockDegrees){double ret=m_rockDegrees;
                                              m_rockDegrees = rockDegrees;
                                              return ret;}
    
    int setRockType(RockType rockType);//{m_rockType = rockType;}
    int setRockType(int rockType);//{m_rockType = rockType;}

    double RAX()const{return m_RAX;}
    double RAZ()const{return m_RAZ;}
    double DECX()const{return m_DECX;}
    double DECZ()const{return m_DECZ;}
    double RAZenith()const{return m_RAZenith;}
    double DECZenith()const{return m_DECZenith;}

    Hep3Vector position(double seconds)const{
        double time = m_earthOrbit->dateFromSeconds(seconds);
        return m_earthOrbit->position(time);
        /*return m_position;*/} //interface to EarthOrbit::position()
    
    protected:
        // singleton - protect ctor/dtor
        GPS();
        virtual ~GPS();
        
        void    time ( GPStime );       // set time
        std::pair<double,double> m_rotangles;  //angles for coordinate rotation (rocking angle)
        
        // friends
        friend class FluxGenerator;
        
    private:
        static GPS* s_instance;
        astro::EarthOrbit* m_earthOrbit; //orbital position object, from the astro package.
         
        double  m_expansion;    // orbit expansion factor
        GPStime m_time;	    // global time
        double  m_sampleintvl;  // interval to sample for each pt. in the orbit - to normalize spectra
        double m_RAX,m_RAZ,m_DECX,m_DECZ; //pointing characteristics.
        double m_RAZenith,m_DECZenith;  //pointing characteristic of the zenith direction.
        Hep3Vector m_position; //current vector position of the LAT.
        // notification
        Subject    m_notification; 
        double m_rockDegrees; //number of degrees to "rock" the spacecraft, along the local x axis.  
        RockType m_rockType;//current rocking scheme
};

inline std::istream&    operator>>(std::istream& i, GPS::Coords& c) {
    double  p, y, r, lt, ln, tm, ph;
    i >> lt; i >> ln; i >> p; i >> y; i >> r; i >> tm; i >> ph;
    c = GPS::Coords( lt, ln, p, y, r, tm, ph );
    return i;
}

inline std::ostream&    operator<<(std::ostream& o, const GPS::Coords& c) {
    o << ' ' << c.lat() << ' ' << c.lon() << ' '  
        << c.time() <<' ' << c.phase();
    return o;
}
#endif // !defined(AFX_GPS_H__F9844433_4E64_11D2_B4DD_00A0C9960210__INCLUDED_)
