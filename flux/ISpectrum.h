#if !defined(AFX_ISPECTRUM_H__9D00078B_7121_473C_813A_827B2F4126A5__INCLUDED_)
#define AFX_ISPECTRUM_H__9D00078B_7121_473C_813A_827B2F4126A5__INCLUDED_

/** 
* \class ISpectrum
*
* \brief The virtual interface for Spectrum-type objects.
*
* Class for holding function definitions of Spectrums...
* an abstract base class.
*
* Notes: The input variable "time" means the
time when the event is called - that is, when the function flux(time1)
gets called, what should be returned is the flux (c/s/m^2/sr) for the
source at that time (time1), and/or when interval(time1) gets called, it
is asking for the interval between time1 and the time when the next
particle arrives.

However, using these finctions is tricky, since they seem to
contradict each other (i.e. "what happens if the flux at a certain time
is really low, but the interval to the next particle is really short?")
The answer lies in understanding which figure the FluxSvc code will use
to find the next particle, and it works this way:

1. FluxSvc will first call interval(time1) to find out when the next particle
comes. if this is a "good" (nonnegative) number, then this will be the
time when the next particle is understood to arrive.  If this is instead
a negative number (say -1), then the process will go to step 2:

2. FluxSvc will then call flux(time1) to determine the average flux for
that time, and then calculate the arrival time for the next particle using
the poisson distribution, as well as the solid angle of the source
(since flux is in units of c/s/m^2/sr.)

So, if interval(time) is set up to already know when the next particle will
arrive (this is necessary for sources with a high time-variance), then
FluxSvc will just use it, and otherwise, flux(time) will be used to
calculate the next time.

* \author Sean Robinson
*
* $Header: /nfs/slac/g/glast/ground/cvs/flux/flux/ISpectrum.h,v 1.2 2004/09/29 18:27:47 cohen Exp $
*/

#include <string>
#include <utility> // for std::pair


class ISpectrum  
{
    
public:
    
    ///  particle name that must be known to the particle service
    virtual const char * particleName()const=0;
    
    /** calculate the flux, particles/m^2/sr. (default zero)
        @param time the mission elapsed time in sec
    */
    virtual double    flux (double time ) const=0;
    
    /// return effective solid angle that will be used to determine the actual rate.
    /// If the source declared in the xml file 
    /// is using the "use_spectrum" tag, telling FluxSvc
    /// to use all the directional information from the class itself, then the
    /// solid_angle function really only tells FluxSvc what angular area the flux
    /// is over, so that it can calculate a rate.
    virtual double solidAngle()const{return  1.0;} //flag that doesen't calculate.
    
    /// return a title describing the spectrum	
    virtual std::string title()const=0;
    
    /*! a (randomized) interval to the next event.  
    @param time the mission elapsed time in sec
     For time-independent rate, should correspond to exponential( 1/rate() )
     Return negative to do this with flux()*solidAngle().
     needs to know the cross-sectional area?
    */
    virtual double interval (double time)=0;
    
    /// return energy, either GeV or MeV
    virtual double energy( double time=0)=0;

    /** return direction in a pair:
    @param energy The generated energy, from previous call to energy
    @return direction is either in the format (cos theta, phi) for
    (zenith-local coordinates, or (l,b) (galactic coordinates).
    */
    virtual std::pair<double,double> dir(double energy)=0;
    
    /** dummy set methods that are actually defaulted in Spectrum class.
     * These are needed to parse info from the XML (see FluxSource class)
     */
    virtual void setParticleName(const std::string& ){;}
    virtual void setInGeV(const bool ){;}
    virtual void setFlux(double ){;}
};


#endif // !defined(AFX_ISPECTRUM_H__9D00078B_7121_473C_813A_827B2F4126A5__INCLUDED_)
