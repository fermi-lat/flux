// GPS.cxx: implementation of the GPS class.
// $Id: GPS.cxx,v 1.19 2003/11/24 18:47:27 srobinsn Exp $
//////////////////////////////////////////////////////////////////////

#include "flux/GPS.h"

//#include "Orbit.h"
#include "CLHEP/Random/RandFlat.h"
#include "astro/EarthCoordinate.h"

#include <iomanip>


// static declaration

GPS*	GPS::s_instance = 0;

GPS::GPS() 
:m_rotangles(std::make_pair<double,double>(0.,0.)),
m_earthOrbit(new astro::EarthOrbit),
m_expansion(1.),    // default expansion:regular orbit for now
m_time(0.), 
m_lastQueriedTime(-1.),
m_sampleintvl(30.), // update position every 30 seconds
m_rockDegrees(35.),
m_rockType(NONE),
m_rockNorth(0)
{}

GPS::Coords::Coords( double alat, double alon, double apitch
                    , double ayaw, double aroll, GPStime atime, double aphase) 
                    : m_time(atime), m_lat(alat), m_lon(alon), 
                    m_pitch(apitch), m_yaw(ayaw), m_roll(aroll), m_phase(aphase)
{}

GPS::Coords::Coords(){}


GPS::~GPS ()
{ }//delete m_orbit; }


void GPS::synch ()
{
    static bool first=true;
    bool changed=  false;
    static GPStime  last_time = time();

    if (Scheduler::instance()->running()) {
        time( Scheduler::instance()->elapsed_time() );
        changed = true; // maybe have threshold before nofitying?
    }

    // If elapsed time exceeds interval then update
    if ((time() - last_time) > m_sampleintvl) {
        last_time = time();
        changed = true;    
    }

    // notify observers if changed (or first time thru)
    if( changed || first) notifyObservers();
    first=false;

}

double GPS::lat()const{return m_lat;}     
/// present longitude    
double GPS::lon()const{return m_lon;}  	
/// pointing characteristics	
double GPS::RAX()const{return m_RAX;}    
double GPS::RAZ()const{return m_RAZ;}    
double GPS::DECX()const{return m_DECX;}    
double GPS::DECZ()const{return m_DECZ;}    
double GPS::RAZenith()const{return m_RAZenith;}    
double GPS::DECZenith()const{return m_DECZenith;}

GPStime	GPS::time ()  const
{ 
    return m_time;
}


double   GPS::expansion () const
{
    return m_expansion;
}

void GPS::pass ( double t )
{ 
    if (!Scheduler::instance()->running())	{
        time(time() + t);
    }   // cannot pass when scheduler is in control!
}

void GPS::expansion ( double e )
{
    m_expansion = e; 
}

void GPS::time ( GPStime t )
{
    m_time = t;
    GPS::instance()->getPointingCharacteristics(t);
}

GPS*	GPS::instance() 
{ return (s_instance != 0) ? s_instance : (s_instance = new GPS()); }

void GPS::kill()
{
    delete s_instance;
    s_instance = 0;
}

void    GPS::sampleintvl ( GPStime s )
{
    m_sampleintvl = s;
}

GPStime  GPS::sampleintvl () const
{
    return m_sampleintvl;
}


void    GPS::printOn(std::ostream& out) const
{
    out << "GPS: time=" << time() 
        << ", lat, lon=" 
        << std::setprecision(4) << lat() << ", " << lon() 
        << std::endl;
}


//access m_rotangles
std::pair<double,double> GPS::rotateAngles(){
    return m_rotangles;

}

//set m_rotangles
void GPS::rotateAngles(std::pair<double,double> coords){
    m_rotangles=coords;
    //m_rockType = EXPLICIT;
}

/// set the desired pointing history file to use:
void GPS::setPointingHistoryFile(std::string fileName){
    m_pointingHistoryFile = fileName;
    m_rockType = HISTORY;
    setUpHistory();
}

HepRotation GPS::transformToGlast(double seconds,CoordSystem index){
    ///Purpose:  return the rotation from a given system to
    ///the satellite system.  this should eventually be the
    ///only function used external to GPS!

    ///Input:  Current time, and a "coordinate system" index:
    ///0: The original direction is in the GLAST frame already
    ///1: rotate from the earth-zenith frame to GLAST
    ///2: rotate from the cartesian celestial coordinate system (like a SkyDir)

    ///Output:  3x3 rocking-angle transformation matrix.
    HepRotation trans;

    if(index==0){
        //do nothing - we are already in the GLAST frame.
    }else if(index==1){
        //earth-zenith to GLAST - just the rocking rotation.
        trans=rockingAngleTransform(seconds);
    }else if(index==2){
        //SkyDir to GLAST - the transformCELToGlast
        trans=transformCelToGlast(seconds);
    }else{
        throw "index out of range in GPS transformtaion matrix function!";
    }
    
    return trans;
}

HepRotation GPS::rockingAngleTransform(double seconds){
    ///Purpose:  return the rotation to correct for satellite rocking.
    ///Input:  Current time
    ///Output:  3x3 rocking-angle transformation matrix.
    using namespace astro;

    // set the needed pointing/location variables:
    getPointingCharacteristics(seconds);

    // now, we want to find the proper transformation for the rocking angles:
    HepRotation rockRot;
    if(m_rockType == EXPLICIT){
        rockRot.rotateX(m_rotangles.first).rotateZ(m_rotangles.second);
    }else{
        //don't do anything - the new pointing characteristics have already been taken account of
        //in getPointingCharacteristics().
        //rockRot.rotateX(m_rockNorth);
        SkyDir dirZenith(m_RAZenith,m_DECZenith,SkyDir::CELESTIAL);
        SkyDir dirXZenith(m_RAXZenith,m_DECXZenith);

    // orthogonalize, since interpolation and transformations destory orthogonality (limit is 10E-8)
    Hep3Vector xhat = dirXZenith() -  (dirXZenith().dot(dirZenith())) * dirZenith() ;
    //so now we know where the x and z axes of the zenith-pointing frame point in the celestial frame.
    //what we want now is to make cel, where
    //cel is the matrix which rotates (cartesian)local coordinates into (cartesian)celestial ones
    HepRotation cel(xhat , dirZenith().cross(xhat) , dirZenith());
    HepRotation temp = transformCelToGlast(seconds) * cel;
    rockRot=temp;
    }

    return rockRot;
}

HepRotation GPS::CELTransform(double seconds){
    /// Purpose:  Return the 3x3 matrix which transforms a vector from a galactic 
    /// coordinate system to a local coordinate system.
    //THIS FUNCTION SHOULD BE REMOVED!  IT IS ONLY USED IN A FEW PLACES!!
    double degsPerRad = 180./M_PI;

    // set the needed pointing/location variables:
    getPointingCharacteristics(seconds);

    //cel is the matrix which rotates (cartesian)local coordinates into (cartesian)celestial ones
    HepRotation cel(transformCelToGlast(seconds).inverse());

    //std::cout << "time is " << seconds << std::endl;
    //m_orbit->displayRotation(cel);

    //gal is the matrix which rotates (cartesian)celestial coordiantes into (cartesian)galactic ones
    HepRotation gal;
    gal.rotateZ(-282.8592/degsPerRad).rotateX(-62.8717/degsPerRad).rotateZ(32.93224/degsPerRad);
    //so gal*cel should be the matrix that makes local coordiates into galactic ones.
    HepRotation glstToGal=gal*cel;
    return glstToGal.inverse();
}

HepRotation GPS::transformCelToGlast(double seconds){
    /// Purpose:  Return the 3x3 matrix which transforms a vector from a celestial 
    /// coordinate system (like a SkyDir vector) to a local coordinate system (like the FluxSvc output).
    using namespace astro;
    //double degsPerRad = 180./M_PI;

    // set the needed pointing/location variables:
    getPointingCharacteristics(seconds);

    double time = m_earthOrbit->dateFromSeconds(seconds);

    //m_position = m_earthOrbit->position(time);

    SkyDir dirZ(m_RAZ,m_DECZ,SkyDir::CELESTIAL);
    SkyDir dirX(m_RAX,m_DECX);

    // orthogonalize, since interpolation and transformations destory orthogonality (limit is 10E-8)
    Hep3Vector xhat = dirX() -  (dirX().dot(dirZ())) * dirZ() ;
    //so now we know where the x and z axes of the zenith-pointing frame point in the celestial frame.
    //what we want now is to make cel, where
    //cel is the matrix which rotates (cartesian)local coordinates into (cartesian)celestial ones
    HepRotation cel(xhat , dirZ().cross(xhat) , dirZ());
    return cel.inverse();
}

HepRotation GPS::transformGlastToGalactic(double seconds){
    return (CELTransform(seconds).inverse())*(rockingAngleTransform(seconds).inverse());
}

void GPS::getPointingCharacteristics(double inputTime){
    //this function calculates all the relevant satellite data. Note that the rocking is done here as well, 
    //so that pointing characteristics come out right.
    using namespace astro;

    //decide if the time has changed;  if it has not, we have already calculated all of the following information:
    if(m_lastQueriedTime==inputTime)return;

    //if not, set the last time to this one:
    m_lastQueriedTime=inputTime;

    //and then get the appropriate julian date:
    double time = m_earthOrbit->dateFromSeconds(inputTime);

    double inclination = m_earthOrbit->inclination();
    double orbitPhase = m_earthOrbit->phase(time);
    m_position = m_earthOrbit->position(time);

    //first make the directions for the x and Z axes, as well as the zenith direction, and latitude/longitude
    double lZ,bZ,raX,decX;
    //before rotation, the z axis points along the zenith:
    if(m_rockType == POINT){
        lZ=m_rotangles.first;
        bZ=m_rotangles.second;
        SkyDir tempDirZ(lZ,bZ,astro::SkyDir::GALACTIC);
        raX = tempDirZ.ra()-90.0;
        decX = 0.;
        astro::EarthCoordinate earthpos(m_position,time);
        m_lat = earthpos.latitude();
        m_lon = earthpos.longitude();
        //now set the zenith direction before the rocking.
        m_RAZenith = tempDirZ.ra();
        m_DECZenith = tempDirZ.dec();
    }else if(m_rockType == HISTORY){
        setInterpPoint(inputTime);
        SkyDir dirZenith(m_currentInterpPoint.position.unit());
        SkyDir dirZ(m_currentInterpPoint.dirZ);
        SkyDir dirX(m_currentInterpPoint.dirX);
        lZ=dirZ.l();
        bZ=dirZ.b();
        raX=dirX.ra();
        decX=dirX.dec();
        m_lat = m_currentInterpPoint.lat;
        m_lon = m_currentInterpPoint.lon;
        //now set the zenith direction before the rocking.
        m_RAZenith = dirZenith.ra();
        m_DECZenith = dirZenith.dec();
    }else{
        //ok, get the pointing from earthOrbit.
        SkyDir tempDirZ(m_position.unit());
        lZ=tempDirZ.l();
        bZ=tempDirZ.b();
        raX = tempDirZ.ra()-90.0;
        decX = 0.;
        astro::EarthCoordinate earthpos(m_position,time);
        m_lat = earthpos.latitude();
        m_lon = earthpos.longitude();
        //now set the zenith direction before the rocking.
        m_RAZenith = tempDirZ.ra();
        m_DECZenith = tempDirZ.dec();
    }

    SkyDir dirZ(lZ,bZ,SkyDir::GALACTIC);
    SkyDir dirX(raX,decX);
    //before rotation, the z axis points along the zenith:
    SkyDir dirZenith(dirZ.dir());
    //also, before Rotation, set the characteristics of the zenith x-direction:
    m_RAXZenith = raX;
    m_DECXZenith = decX;

    if(m_rockType != HISTORY){
        //rotate the x direction so that the x direction points along the orbital direction.
        dirX().rotate(dirZ.dir() , inclination*cos(orbitPhase));
    }

    // now, we want to find the proper transformation for the rocking angles:
    //HepRotation rockRot(Hep3Vector(0,0,1).cross(dirZ.dir()) , m_rockNorth);    
    //and apply the transformation to dirZ and dirX:
    m_rockNorth = m_rockDegrees*M_PI/180;
    //here's where we decide how much to rock about the x axis.  this depends on the 
    //rocking mode.
    if(m_rockType == NONE){
        m_rockNorth = 0.;
    }else if(m_rockType == UPDOWN){
        if(m_DECZenith <= 0) m_rockNorth *= -1.;
    }else if(m_rockType == SLEWING){
        //slewing is experimental
        if(m_DECZenith <= 0) m_rockNorth *= -1.;
        if(m_DECZenith >= -5.5 && m_DECZenith <= 5.5){
            m_rockNorth -= m_rockNorth*((5.5-fabs(m_DECZenith))/5.5);
        }
    }else if(m_rockType == ONEPERORBIT){
        while(orbitPhase >2.*M_2PI) {orbitPhase -= 2.*M_2PI;}
        if(orbitPhase <= M_2PI) m_rockNorth *= -1.;
    }else{
        //important - this includes EXPLICIT rocking angles - they
        //are currently still handled by rockingTransform().
        m_rockNorth = 0.;
    }

    dirZ().rotate(dirX.dir() , m_rockNorth);

    m_RAX = dirX.ra();
    m_RAZ = dirZ.ra();
    m_DECX = dirX.dec();
    m_DECZ = dirZ.dec();

    //a test - to ensure the rotation went properly
    //std::cout << " degrees between xhat and zhat directions: " <<
    //    dirZ.difference(dirX)*180./M_PI << std::endl;

}

int GPS::setRockType(int rockType){
    //get the current state
    int ret;
    if(m_rockType == NONE)ret = 0;
    if(m_rockType == UPDOWN)ret = 1;
    if(m_rockType == SLEWING)ret = 2;
    if(m_rockType == ONEPERORBIT)ret = 3;
    if(m_rockType == EXPLICIT)ret = 4;
    if(m_rockType == POINT)ret = 5;
    if(m_rockType == HISTORY)ret = 6;


    m_rockType = NONE;
    if(rockType == 1) m_rockType = UPDOWN;
    if(rockType == 2) m_rockType = SLEWING;
    if(rockType == 3) m_rockType = ONEPERORBIT;
    if(rockType == 4) m_rockType = EXPLICIT;
    if(rockType == 5) m_rockType = POINT;
    if(rockType == 6) m_rockType = HISTORY;

    return ret;
}

int GPS::setRockType(RockType rockType){
    //get the current state
    int ret;
    if(m_rockType == NONE)ret = 0;
    if(m_rockType == UPDOWN)ret = 1;
    if(m_rockType == SLEWING)ret = 2;
    if(m_rockType == ONEPERORBIT)ret = 3;
    if(m_rockType == EXPLICIT)ret = 4;
    if(m_rockType == POINT)ret = 5;
    if(m_rockType == HISTORY)ret = 6;

    m_rockType = rockType;
    return ret;
}

void GPS::setUpHistory(){
    std::ifstream input_file;
    input_file.open(m_pointingHistoryFile.c_str());

    if(false == input_file.is_open())
    {
        std::cerr << "ERROR:  Unable to open:  " << m_pointingHistoryFile.c_str() << std::endl;
        exit(0);
    }
    else
    {
        double intrvalstart,posx,posy,posz,raz,decz,rax,decx,razenith,deczenith,lon,lat,alt;
        //initialize the key structure:
        while (!input_file.eof()){
            input_file >> intrvalstart;
            input_file >> posx;
            input_file >> posy;
            input_file >> posz;
            input_file >> raz;
            input_file >> decz;
            input_file >> rax;
            input_file >> decx;
            input_file >> razenith;
            input_file >> deczenith;
            input_file >> lon;
            input_file >> lat;
            input_file >> alt;

            POINTINFO temp;
            temp.dirZ=astro::SkyDir(raz,decz);
            temp.dirX=astro::SkyDir(rax,decx);
            temp.lat=lat;
            temp.lon=lon;
            temp.position=Hep3Vector(posx,posy,posz);

            m_pointingHistory[intrvalstart]=temp;

        }
    }
}

void GPS::setInterpPoint(double time){
    //is the current time out of range?
    bool timeTooEarly=false;
    bool timeTooLate=false;
    std::map<double,POINTINFO>::const_iterator iter=m_pointingHistory.upper_bound(time);
    if((time< (*(m_pointingHistory.begin())).first )){
        timeTooEarly=true;
        std::cerr << "WARNING: Time (" << time << ") out of range of times in the pointing database - closest record being used." << std::endl;
        throw "Time out of Range!";
    }else if(iter==m_pointingHistory.end()){
        timeTooLate=true;
        std::cerr << "WARNING: Time (" << time << ") out of range of times in the pointing database - closest record being used." << std::endl;
        throw "Time out of Range!";
    }
    //get the point after "time"
    double rax2=(*iter).second.dirX.ra();
    double decx2=(*iter).second.dirX.dec();
    double raz2=(*iter).second.dirZ.ra();
    double decz2=(*iter).second.dirZ.dec();
    double lat2=(*iter).second.lat;
    double lon2=(*iter).second.lon;
    Hep3Vector pos2=(*iter).second.position;
    double time2=(*iter).first;
	astro::SkyDir dirZ2=(*iter).second.dirZ;
	astro::SkyDir dirX2=(*iter).second.dirX;

    //then get the details from the previous point:
    iter--;
    double rax1=(*iter).second.dirX.ra();
    double decx1=(*iter).second.dirX.dec();
    double raz1=(*iter).second.dirZ.ra();
    double decz1=(*iter).second.dirZ.dec();
    double lat1=(*iter).second.lat;
    double lon1=(*iter).second.lon;
    Hep3Vector pos1=(*iter).second.position;
    double time1=(*iter).first;
	astro::SkyDir dirZ1=(*iter).second.dirZ;
	astro::SkyDir dirX1=(*iter).second.dirX;

    //the proportional distance between the first point and the interpolated point
    double prop;
    if(timeTooEarly){
        //use the later point(first in database)
        prop=1.0;
    }else if(timeTooLate){
        //use the earlier point(last in database)
        prop=0.0;
    }else{
        //do the interpolation
        prop= 1.0 - ((time2-time)/(time2-time1));
    }

    m_currentInterpPoint.position=pos1+((pos2-pos1)*prop);

    m_currentInterpPoint.lat=lat1+((lat2-lat1)*prop);
    m_currentInterpPoint.lon=lon1+((lon2-lon1)*prop);	
    m_currentInterpPoint.dirZ=astro::SkyDir(raz1+((raz2-raz1)*prop),decz1+((decz2-decz1)*prop));
    m_currentInterpPoint.dirX=astro::SkyDir(rax1+((rax2-rax1)*prop),decx1+((decx2-decx1)*prop));

	//this piece of code should just handle the "wraparound" cases:
	if(fabs(lon1-lon2) >= 330.){
		//we have gone off one end of the longitude scale.
		double lonlesser=std::max(lon1,lon2);
		double longreater=std::min(lon1,lon2)+360.;
		double interpLon = lonlesser+((longreater-lonlesser)*prop);
		while(interpLon > 360.)interpLon -= 360.;
		m_currentInterpPoint.lon=interpLon;
	}

	//now, find dirZ and X between the two nearest cases.
	m_currentInterpPoint.dirZ=astro::SkyDir(dirZ1()+((dirZ2()-dirZ1())*prop));
    m_currentInterpPoint.dirX=astro::SkyDir(dirX1()+((dirX2()-dirX1())*prop));

    //now regenerate X to be perpindicular to Z (ir should be almost there anyway).
    Hep3Vector dirY( m_currentInterpPoint.dirZ().cross(m_currentInterpPoint.dirX()) );
    m_currentInterpPoint.dirX=astro::SkyDir( dirY.cross(m_currentInterpPoint.dirZ()) );

}