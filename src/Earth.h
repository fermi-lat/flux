/**
 * @file Earth.h
 * @brief A source class for the flux package that 
 * simulates the gamma emission of the Earth based on
 * measurements made by EGRET 
 * (see http://xxx.lanl.gov/abs/astro-ph/0410487 )
 * @author D. Petry
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/flux/src/Earth.h,v 1.1 2004/12/21 03:46:41 burnett Exp $
 */

#ifndef mySpectrum_MapSource_h
#define mySpectrum_MapSource_h

#include <vector>

#include "flux/Spectrum.h"

/**
 * @class Earth
 *
 * @brief 
 * A source class for the flux package that
 * simulates the gamma emission of the Earth based on
 * measurements made by EGRET
 * (see http://xxx.lanl.gov/abs/astro-ph/0410487 )
 * @author D. Petry
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/flux/src/Earth.h,v 1.1 2004/12/21 03:46:41 burnett Exp $
 */

class Earth : public Spectrum {

public:

   /// This constructor is required and used in FluxSource for
   /// "SpectrumClass" sources.
   Earth(const std::string &params);

   virtual ~Earth(){}

 
   /// @return Particle type, "gamma".
   virtual const char * particleName() const {return "gamma";}

   /// @return Total flux (photons/m^2).
   /// @param time Simulation time in seconds.
   virtual double flux(double time) const;

   /// @return "Effective" solid angle (sr).
   virtual double solidAngle() const;

   /// @return Title describing the spectrum.
   virtual std::string title() const {return "Earth";}

   /// @return Interval to the next event (seconds)
   virtual double interval(double time);

   /// @return Photon energy (MeV).
   virtual double energy(double time);

   /// @return Photon direction in (zenith angle (degrees), azimuth (degrees,East=0,North=90deg)).
   virtual std::pair<double, double> dir(double energy);

   /// @return energy and photon direction (zenith angle (degrees), azimuth (degrees,East=0,North=90deg)).
   std::pair<double, std::pair<double, double> > photon();

   /// initialize the model parameters assuming altitude "alt" (km)
   void init(double alt, double emin, double emax);

   /// Integral earth flux over all sky from emin to emax (set in init) using isteps integration steps
   double ff(int isteps) const; 

   /// Integral earth flux (in cm^-2s^-1) over part of the sky from emin to emax using isteps integration steps
   double ffpartial(int isteps, 
		    double zamin, double zamax, 
		    double azmin, double azmax) const;

private:

   double m_alt, m_emin, m_emax, m_a[28], m_ftot; // model parameters, set in init()
   int m_version; // the version of the model parameters 
                  //  (used to check if the total flux can be taken from a precomputed value) 
   bool m_eCalled;
   double m_e, m_t, m_p;

   double tp(double ee) const; // Peak position (degrees)
   double sigma(double ee) const; // Peak Width (degrees)
                            //corrected for EGRET PSF, orbital decay, orbital interpolation
   double sigmaUncorrected(double ee) const; // Peak Width as measured by EGRET
   double t0(double ee) const; // Transition point from Gaussian to Exponential 
   double g0(double e) const; // Spectrum of the constant component independent of AZ
   double g1(double e) const; // Spectrum of the component dependent on AZ
   double sigmaz(double t) const; // Std Dev of the AZ dependent component
   double f0(double t, double p, double e) const; // the part zt ZA < t0
   double normcorr(double e) const; // correction of the normalization of fa(t,p,e)
   double fa(double t, double p, double e) const;
   double b(double e) const; // central flux
   double c(double e, double p) const;
   double c0(double e, double p) const;
   double fb(double t, double p, double e) const;
   double f(double t, double p, double e) const;

   double pp(double t, double p, double e) const; // the probability density normalized to give 1 when 
                                            // integrated over the whole sky
   double qn(double e) const; // the normalization of the envelope function
   double q(double e) const; // the envelope function
   double qq(double e) const; // the integral of the envelope function
   double qqinv(double x) const; // the inverse of the integral of q
   void earth(double &t, double &p, double &e) const ; // return zenith angle (deg), azimuth (deg) and energy (MeV)


};

#endif 

