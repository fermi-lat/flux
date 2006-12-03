/** @file SourceDirection.h
@brief SourceDirection declaration

$Header$

*/
#ifndef flux_SourceDirection_h
#define flux_SourceDirection_h

#include "flux/LaunchDirection.h"

#include <string>

class ISpectrum;

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
/** @class SourceDirection
@brief launch strategy derived class

Gets a direction from the ISpectrum class
*/
class SourceDirection : public LaunchDirection{ 
public:
    /** Ctor:
    @param spectrum pointer to the ISpectrum object that will provide the direction
    @param galactic if true, interpret pair as l,b (in degrees); otherwise costh, phi
    */
    SourceDirection(ISpectrum* spectrum, std::string frame );

    void execute(double ke, double time);

    //! solid angle
    virtual double solidAngle()const;

    virtual std::string title()const { return "(use_spectrum frame=\""+m_frameName+"\")";
    }

    /// return the cosine of the angle between the incoming direction and the earth's zenith
    virtual double zenithCosine()const{return m_zenithCos;}

private:
    typedef enum  { ZENITH, EQUATORIAL, GALACTIC, SUN, MOON, INVALID } Frame;
    Frame m_frame;
    ISpectrum* m_spectrum;
    double m_zenithCos;
    std::string m_frameName;

    CLHEP::Hep3Vector solarSystemDir(Frame f, double time);

};

#endif

