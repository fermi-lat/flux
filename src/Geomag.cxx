#include "Geomag.h"
#include <iostream>

#include "Geomag.inc"

double Geomag::L(double lat, double lon) {
    return Geomag::geoInterp(lat, lon, Lvals);
}

double Geomag::B(double lat, double lon) {
    return Geomag::geoInterp(lat, lon, Bvals);
}

double Geomag::geolat(double lat, double lon) {
    return Geomag::geoInterp(lat, lon, glats);
}

double Geomag::geolon(double lat, double lon) {
    return Geomag::geoInterp(lat, lon, glons);
}

double Geomag::L(std::pair<double,double>coords) {
    return Geomag::L(coords.first, coords.second);
}

double Geomag::B(std::pair<double,double>coords) {
    return Geomag::B(coords.first, coords.second);
}

double Geomag::geolat(std::pair<double,double>coords) {
    return Geomag::geolat(coords.first, coords.second);
}

double Geomag::geolon(std::pair<double,double>coords) {
    return Geomag::geolon(coords.first, coords.second);
}

// Allowed Latitude range is -30 to +30
// Allowed Longitude range is -180 to 180 (but 0 to 360 works too)
double Geomag::geoInterp(double lat, double lon, const double * array) {
   int ilat, ilon;
   double a, b;

   if(lon < 0)
      lon += 360;

    ilat = static_cast<int>(lat/5.+6);
    ilon = static_cast<int>(lon/5.);
    a = fmod(lat+30., 5.)/5.;
    b = fmod(lon, 5.)/5.;

	if(abs(lat) > 30)
	{
		std::cout << "Warning -- Geomag::geoInterp -- Latitude is out of range."  << std::endl;
	}

    return array[ilat   + 13*ilon    ] * (1.-a) * (1.-b) +
        array[ilat   + 13*(ilon+1)] * (1.-a) * b      +
        array[ilat+1 + 13*ilon    ] * a      * (1.-b) +
        array[ilat+1 + 13*(ilon+1)] * a      * b      ;
}
