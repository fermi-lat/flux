// GPS.cxx: implementation of the GPS class.
// $Id: GPS.cxx,v 1.7 2003/08/29 09:08:16 srobinsn Exp $
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
m_sampleintvl(30.), // update position every 30 seconds
m_rockDegrees(35.),
m_rockType(NONE)
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

double GPS::lat () const
{
	//before anything, check to see if we are using a history file:
	if(m_rockType == HISTORY){
		instance()->setInterpPoint(time());
		return m_currentInterpPoint.lat;
	}

	double currentTime = time(); 
	double julianDate = m_earthOrbit->dateFromSeconds(currentTime);

	//and here the pointing characteristics of the LAT.
	GPS::instance()->getPointingCharacteristics(currentTime);
	Hep3Vector location = position(currentTime);

	astro::EarthCoordinate earthpos(location,julianDate);
	return earthpos.latitude();
}

double	GPS::lon () const
{ 
	//before anything, check to see if we are using a history file:
	if(m_rockType == HISTORY){
	instance()->setInterpPoint(time());
	return m_currentInterpPoint.lon;
	}

	double currentTime = time();
	double julianDate = m_earthOrbit->dateFromSeconds(currentTime);

	//and here the pointing characteristics of the LAT.
	GPS::instance()->getPointingCharacteristics(currentTime);
	Hep3Vector location = position(currentTime);

	astro::EarthCoordinate earthpos(location,julianDate);
	return earthpos.longitude();
}

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


HepRotation GPS::rockingAngleTransform(double seconds){
	///Purpose:  return the rotation to correct for satellite rocking.
	///Input:  Current time
	///Output:  3x3 rocking-angle transformation matrix.
	using namespace astro;

	double time = m_earthOrbit->dateFromSeconds(seconds);

	double inclination = m_earthOrbit->inclination();
	double orbitPhase = m_earthOrbit->phase(time);
	m_position = m_earthOrbit->position(time);

	double lZ,bZ,raX,decX;
	//before anything, check to see if we are using a history file:
	if(m_rockType == HISTORY){
		setInterpPoint(seconds);
		lZ=m_currentInterpPoint.dirZ.l();
		bZ=m_currentInterpPoint.dirZ.b();
		raX=m_currentInterpPoint.dirX.ra();
		decX=m_currentInterpPoint.dirX.dec();
	}else{
	SkyDir tempdirZ(m_position.unit());
	lZ = tempdirZ.l();
	bZ = tempdirZ.b();
	raX = tempdirZ.ra()-90.;
	decX=0.0;
	}

	SkyDir dirZ(lZ,bZ,SkyDir::GALACTIC);
	SkyDir dirX(raX,decX);

	//rotate the x direction so that the x direction points along the orbital direction.
	dirX().rotate(dirZ.dir() , inclination*cos(orbitPhase));

	//now set the zenith direction to do the rocking properly.
	m_RAZenith = dirZ.ra();
	m_DECZenith = dirZ.dec();

	double rockNorth = m_rockDegrees*M_PI/180;

	//here's where we decide how much to rock about the x axis.  this depends on the 
	//rocking mode.
	if(m_rockType == NONE){
		rockNorth = 0.;
	}else if(m_rockType == UPDOWN){
		if(m_DECZenith <= 0) rockNorth *= -1.;
	}else if(m_rockType == SLEWING){
		//slewing is experimental
		if(m_DECZenith <= 0) rockNorth *= -1.;
		if(m_DECZenith >= -5.5 && m_DECZenith <= 5.5){
			rockNorth -= rockNorth*((5.5-fabs(m_DECZenith))/5.5);
		}
	}else if(m_rockType == ONEPERORBIT){
		while(orbitPhase >2.*M_2PI){ orbitPhase -= 2.*M_2PI;}
		if(orbitPhase <= M_2PI) rockNorth *= -1.;
	}else{
		//for safety and EXPLICIT and POINT
		rockNorth = 0.;
	}
	// now, we want to find the proper transformation for the rocking angles:
	HepRotation rockRot;
	if(m_rockType == EXPLICIT){
		rockRot.rotateX(m_rotangles.first).rotateZ(m_rotangles.second);
	}else{
		//just rock north or south by the desired amount.
		rockRot./*rotateZ(inclination*cos(orbitPhase)).*/rotateX(rockNorth);
	}

	return rockRot;
}

HepRotation GPS::CELTransform(double seconds){
	/// Purpose:  Return the 3x3 matrix which transforms a vector from a galactic 
	/// coordinate system to a local coordinate system.
	using namespace astro;
	double degsPerRad = 180./M_PI;
	HepRotation gal;//,cel;
	double time = m_earthOrbit->dateFromSeconds(seconds);

	m_position = m_earthOrbit->position(time);

	//first make the directions for the x and Z axes, as well as the zenith direction.
	double lZ,bZ,raX,decX;
	//before rotation, the z axis points along the zenith:
	if(m_rockType == POINT){
		lZ=m_rotangles.first;
		bZ=m_rotangles.second;
		SkyDir tempDirZ(lZ,bZ,astro::SkyDir::GALACTIC);
		raX = tempDirZ.ra()-90.0;
		decX = 0.;
	}else if(m_rockType == HISTORY){
		setInterpPoint(seconds);
		SkyDir dirZ(m_currentInterpPoint.dirZ);
		SkyDir dirX(m_currentInterpPoint.dirX);
		lZ=dirZ.l();
		bZ=dirZ.b();
		raX=dirX.ra();
		decX=dirX.dec();
	}else{
		//ok, get the pointing from earthOrbit.
		SkyDir tempDirZ(m_position.unit());
		lZ=tempDirZ.l();
		bZ=tempDirZ.b();
		raX = tempDirZ.ra()-90.0;
		decX = 0.;
	}

	SkyDir dirZ(lZ,bZ,SkyDir::GALACTIC);
	SkyDir dirX(raX,decX);

	//so now we know where the x and z axes of the zenith-pointing frame point in the celestial frame.
	//what we want now is to make cel, where
	//cel is the matrix which rotates (cartesian)local coordinates into (cartesian)celestial ones
	HepRotation cel(dirX() , dirZ().cross(dirX()) , dirZ());

	//std::cout << "time is " << seconds << std::endl;
	//m_orbit->displayRotation(cel);

	//gal is the matrix which rotates (cartesian)celestial coordiantes into (cartesian)galactic ones
	gal.rotateZ(-282.25/degsPerRad).rotateX(-62.6/degsPerRad).rotateZ(33./degsPerRad);
	//so gal*cel should be the matrix that makes local coordiates into galactic ones.
	HepRotation glstToGal=gal*cel;
	return glstToGal.inverse();

}

HepRotation GPS::transformCelToGlast(double seconds){
	/// Purpose:  Return the 3x3 matrix which transforms a vector from a celestial 
	/// coordinate system (like a SkyDir vector) to a local coordinate system (like the FluxSvc output).
	using namespace astro;
	double degsPerRad = 180./M_PI;

	double time = m_earthOrbit->dateFromSeconds(seconds);

	m_position = m_earthOrbit->position(time);

	//first make the directions for the x and Z axes, as well as the zenith direction.
	double lZ,bZ,raX,decX;
	//before rotation, the z axis points along the zenith:
	if(m_rockType == POINT){
		lZ=m_rotangles.first;
		bZ=m_rotangles.second;
		SkyDir tempDirZ(lZ,bZ,astro::SkyDir::GALACTIC);
		raX = tempDirZ.ra()-90.0;
		decX = 0.;
	}else if(m_rockType == HISTORY){
		setInterpPoint(seconds);
		SkyDir dirZ(m_currentInterpPoint.dirZ);
		SkyDir dirX(m_currentInterpPoint.dirX);
		lZ=dirZ.l();
		bZ=dirZ.b();
		raX=dirX.ra();
		decX=dirX.dec();
	}else{
		//ok, get the pointing from earthOrbit.
		SkyDir tempDirZ(m_position.unit());
		lZ=tempDirZ.l();
		bZ=tempDirZ.b();
		raX = tempDirZ.ra()-90.0;
		decX = 0.;
	}

	SkyDir dirZ(lZ,bZ,SkyDir::GALACTIC);
	SkyDir dirX(raX,decX);

	//so now we know where the x and z axes of the zenith-pointing frame point in the celestial frame.
	//what we want now is to make cel, where
	//cel is the matrix which rotates (cartesian)local coordinates into (cartesian)celestial ones
	HepRotation cel(dirX() , dirZ().cross(dirX()) , dirZ());
	return cel.inverse();
}

HepRotation GPS::transformGlastToGalactic(double seconds){
	return (CELTransform(seconds).inverse())*(rockingAngleTransform(seconds).inverse());
}

void GPS::getPointingCharacteristics(double seconds){
	//this is being used by exposureAlg right now, and should be reworked
	//to use the rest of this class.
	using namespace astro;

	double time = m_earthOrbit->dateFromSeconds(seconds);

	double inclination = m_earthOrbit->inclination();
	double orbitPhase = m_earthOrbit->phase(time);
	m_position = m_earthOrbit->position(time);

	//first make the directions for the x and Z axes, as well as the zenith direction.
	double lZ,bZ,raX,decX;
	//before rotation, the z axis points along the zenith:
	if(m_rockType == POINT){
		lZ=m_rotangles.first;
		bZ=m_rotangles.second;
		SkyDir tempDirZ(lZ,bZ,astro::SkyDir::GALACTIC);
		raX = tempDirZ.ra()-90.0;
		decX = 0.;
	}else if(m_rockType == HISTORY){
		setInterpPoint(seconds);
		SkyDir dirZ(m_currentInterpPoint.dirZ);
		SkyDir dirX(m_currentInterpPoint.dirX);
		lZ=dirZ.l();
		bZ=dirZ.b();
		raX=dirX.ra();
		decX=dirX.dec();
	}else{
		//ok, get the pointing from earthOrbit.
		SkyDir tempDirZ(m_position.unit());
		lZ=tempDirZ.l();
		bZ=tempDirZ.b();
		raX = tempDirZ.ra()-90.0;
		decX = 0.;
	}

	SkyDir dirZ(lZ,bZ,SkyDir::GALACTIC);
	SkyDir dirX(raX,decX);
	//before rotation, the z axis points along the zenith:
	SkyDir dirZenith(dirZ.dir());

	if(m_rockType != HISTORY){
	//rotate the x direction so that the x direction points along the orbital direction.
	dirX().rotate(dirZ.dir() , inclination*cos(orbitPhase));
	}

	//now set the zenith direction before the rocking.
	m_RAZenith = dirZ.ra();
	m_DECZenith = dirZ.dec();

	// now, we want to find the proper transformation for the rocking angles:
	//HepRotation rockRot(Hep3Vector(0,0,1).cross(dirZ.dir()) , rockNorth);    
	//and apply the transformation to dirZ and dirX:
	double rockNorth = m_rockDegrees*M_PI/180;
	//here's where we decide how much to rock about the x axis.  this depends on the 
	//rocking mode.
	if(m_rockType == NONE){
		rockNorth = 0.;
	}else if(m_rockType == UPDOWN){
		if(m_DECZenith <= 0) rockNorth *= -1.;
	}else if(m_rockType == SLEWING){
		//slewing is experimental
		if(m_DECZenith <= 0) rockNorth *= -1.;
		if(m_DECZenith >= -5.5 && m_DECZenith <= 5.5){
			rockNorth -= rockNorth*((5.5-fabs(m_DECZenith))/5.5);
		}
	}else if(m_rockType == ONEPERORBIT){
		while(orbitPhase >2.*M_2PI) {orbitPhase -= 2.*M_2PI;}
		if(orbitPhase <= M_2PI) rockNorth *= -1.;
	}else{
		rockNorth = 0.;
	}

	dirZ().rotate(dirX.dir() , rockNorth);

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
	std::map<double,POINTINFO>::const_iterator iter=m_pointingHistory.upper_bound(time);
	if((time< (*(m_pointingHistory.begin())).first )||iter==m_pointingHistory.end()) std::cerr << "WARNING: Time out of scope of pointing database" << std::endl;
	//get the point after "time"
	double rax2=(*iter).second.dirX.ra();
	double decx2=(*iter).second.dirX.dec();
	double raz2=(*iter).second.dirZ.ra();
	double decz2=(*iter).second.dirZ.dec();
	double lat2=(*iter).second.lat;
	double lon2=(*iter).second.lon;
	Hep3Vector pos2=(*iter).second.position;
	double time2=(*iter).first;	
	
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

	//the proportional distance between the first point and the interpolated point
	double prop=1.0 - ((time2-time)/(time2-time1));

	m_currentInterpPoint.position=pos1+((pos2-pos1)*prop);
	m_currentInterpPoint.lat=lat1+((lat2-lat1)*prop);
	m_currentInterpPoint.lon=lon1+((lon2-lon1)*prop);	
	m_currentInterpPoint.dirZ=astro::SkyDir(raz1+((raz2-raz1)*prop),decz1+((decz2-decz1)*prop));
	m_currentInterpPoint.dirX=astro::SkyDir(rax1+((rax2-rax1)*prop),decx1+((decx2-decx1)*prop));
}