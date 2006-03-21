/**
 * @file Earth.cxx
 * @brief A phenomenological model of the Earth based on EGRET measurements
 * @author D. Petry
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/flux/src/Earth.cxx,v 1.5 2005/05/26 03:53:33 burnett Exp $
 */

#include <iostream>

#include <cmath>
#include <cstdlib>
#include <stdexcept>

#include "facilities/Util.h"

#include "flux/SpectrumFactory.h"
#include "flux/EventSource.h"
#include "CLHEP/Random/RandFlat.h"


#include "Earth.h"

static SpectrumFactory<Earth> factory;
const ISpectrumFactory& EarthFactory = factory;


double Earth::tp() const {
// Peak position
    return m_a[1];
}

double Earth::sigma() const {
// Peak Width (deg) corrected for EGRET PSF, orbital decay, orbital interpolation
    return m_a[3]/m_a[4]/4.*180./M_PI * m_a[2];
}

double Earth::t0() const {
// Transition point from Gaussian to Exponential 
    return tp() + 1.5 * sigma() * m_a[2];  
}

double Earth::g0(double e) const {
// Spectrum of the constant component independent of AZ
    return m_a[5]*pow(e,m_a[6])*exp(-e/m_a[7]);
}

double Earth::g1(double e) const {
// Spectrum of the component dependent on AZ
    return m_a[8]*pow(e,m_a[9])*exp(-e/m_a[10]);
}

double Earth::sigmaz(double t) const {
// Std Dev of the AZ dependent component
    return m_a[11] * exp(-0.5*pow(((t-m_a[12])/m_a[13]),2.));
}

double Earth::f0(double t, double p, double e) const {
// the part zt ZA < t0
    return g0(e) + g1(e) * exp(-0.5*pow(((p-m_a[14])/sigmaz(t)),2.));
}


double Earth::normcorr(double e) const {
// correction of the normalization of fa(t,p,e)
// necessary because we are replacing the measured peak width 
// by sigma() which is corrected for the EGRET PSF and effects
// from orbital decay and orbital interpolation
  double rval, loge;
  loge = log10(e); // e (MeV)
  rval = 1./(0.61881 - 0.16826 * loge + 0.35618E-01 * loge*loge)
      / (0.71811 + 0.48235E-01 * loge) ; 
  if (rval < 1.) {
    rval = 1.;
  }
  return rval; 
}


double Earth::fa(double t, double p, double e) const {
    if(t <= t0()+1E-10  &&  t > t0()-5.*sigma()){ // the addition of 1E-10 is necessary
                                                  // to overcome a strange effect in the
                                                  // g++  3.2.3 optimizer
        return normcorr(e)*exp(-0.5*pow(((t-tp())/sigma()),2.))*f0(t,p,e);
    } else{
        return 0.;
    }
}

double Earth::b(double e) const {
// central flux
    return m_a[15]*pow(e,m_a[16])*exp(-e/m_a[17]);
}

double Earth::c(double e, double p) const {
    return (log(fa(t0(),p,e))-log(b(e)))/(t0()-180.);
}

double Earth::c0(double e, double p) const {
    return log(fa(t0(),p,e))-c(e,p)*t0();
}

double Earth::fb(double t, double p, double e) const {
    if(t > t0()+1E-10){
        return exp(c0(e,p)+c(e,p)*t);
    } 
    else{
        return 0.;
    }
}

double Earth::f(double t, double p, double e) const {
    if(t < 90. || t > 180.){
        return 0.;
    } 
    else{
        return fa(t,p,e) + fb(t,p,e);
    }
}

double Earth::ff(int isteps) const {
// Integral of earth flux over all sky from m_a[21] MeV to m_a[22] (MeV)
//  using isteps integration steps
    int i, j, k;
    double theint, t, p, dt, dp, dloge, domega, dtheint;
    double radpdg;
    double steps, logemin, logemax;
    double f1, f2, e1, e2, t1, t2;
    double norm = 0.;
    double alpha;
// check if precomputed value available
    theint = 0.;
    if(m_version == 18052005 
       && m_a[21] == 20. 
       && m_a[22] >= 55000. 
       && isteps >= 200 
       && m_a[19] == 565.){
	std::cout << "Earth: Using precomputed value for integral flux." << std::endl;
	theint = 0.0833719L;
    }
    else if(m_version == 18052005 
       && m_a[21] == 10. && m_a[22] >= 55000. 
       && isteps >= 200 
       && m_a[19] == 565.){
	std::cout << "Earth: Using precomputed value for integral flux." << std::endl;
	theint = 0.18104L;
    }

    if(theint == 0.){
	std::cout << "Earth: computing value for integral flux ..." << std::endl;

        radpdg = M_PI/180.;  
        steps = isteps/1.;
        logemin = log(m_a[21]);
        logemax = log(m_a[22]);
        dp = 360./steps;
        dloge =(logemax-logemin)/steps;
        dt = 180./steps;

        for(i=0;i<isteps;i++){
            
            if(i/100. == i/100){
                std::cout << i <<" of " << isteps << std::endl;
            }

            e1 = exp(logemin + i*dloge);
            e2 = exp(logemin + (i+1)*dloge);
            
            f1 = f2 = 0.;

            for(j=0;j<isteps;j++){
                t1 = j*dt;
                t2 = (j+1)*dt;
                t = (j+0.5)*dt;
                domega = 2*M_PI*(cos(t1*radpdg)-cos(t2*radpdg))/steps;

                for(k=0;k<isteps;k++){
                    p = (k+0.5)*dp;
                    f1 = f1 + f(t,p,e1)*domega;
                    f2 = f2 + f(t,p,e2)*domega;
                }
            }

            alpha = -(log(f1)-log(f2))/(log(e1)-log(e2));
            if(pow(e1,alpha) == HUGE_VAL) {
                dtheint = 0.;
            }
            else {
                norm = f1/pow(e1,-alpha);
                dtheint = norm * ( 
                    pow(e2,-alpha+1.)/(-alpha+1.)
                    -  pow(e1,-alpha+1.)/(-alpha+1.)
                    );
                theint = theint + dtheint;
            }
            std::cout << "in ff: " 
                      << theint << " "
                      << dtheint << " "
                      << alpha << " "
                      << norm << " "
                      << e1 << " "
                      << e2 << " "
                      << std::endl;
        }
    }
    //    std::cout << "in ff: ftot = " << theint << std::endl; 
    return theint;
}

double Earth::ffpartial(int isteps, 
		 double zamin, double zamax, 
		 double azmin, double azmax) const {
// Integral of earth flux (in cm^-2s^-1) over part of the sky from emin to emax
//  using isteps integration steps
    int i, j, k;
    double t, p, dt, dp, dloge, omega, domega, theint, dtheint;
    double radpdg;
    double steps, logemin, logemax;
    double f1, f2, e1, e2, t1, t2;
    double norm; 
    double alpha;
    double ringFraction;
// check if precomputed value available (not yet implemented)
    theint = 0.;
    omega = 0.;
    if(theint == 0.){
        radpdg = M_PI/180.;  
        steps = isteps/1.;
        logemin = log(m_a[21]);
        logemax = log(m_a[22]);
        dp = (azmax-azmin)/steps;
	ringFraction = dp/360.;
        dloge =(logemax-logemin)/steps;
        dt = (zamax-zamin)/steps;

        for(i=0;i<isteps;i++){ // loop over energy
            
            if(i/100. == i/100){
                std::cout << i <<" of " << isteps << std::endl;
            }

            e1 = exp(logemin + i*dloge);
            e2 = exp(logemin + (i+1)*dloge);
            
            f1 = f2 = 0.;

            for(j=0;j<isteps;j++){ // loop over zenith angle
                t1 = j*dt+zamin;
                t2 = (j+1)*dt+zamin;
                t = (j+0.5)*dt+zamin;
                domega = 2*M_PI*(cos(t1*radpdg)-cos(t2*radpdg))*ringFraction;
                for(k=0;k<isteps;k++){ // loop over azimuth
                    p = (k+0.5)*dp + azmin;
                    f1 = f1 + f(t,p,e1)*domega;
                    f2 = f2 + f(t,p,e2)*domega;
                    omega = omega + domega;
                }
            }

            alpha = -(log(f1)-log(f2))/(log(e1)-log(e2));
            if(pow(e1,alpha) == HUGE_VAL) {
                dtheint = 0.;
            }
            else {
                norm = f1/pow(e1,-alpha);
                dtheint = norm * ( 
                    pow(e2,-alpha+1.)/(-alpha+1.)
                    -  pow(e1,-alpha+1.)/(-alpha+1.)
                    );
                theint = theint + dtheint;
            }
        }
    }
    omega = omega/steps;
    return theint;
}


void Earth::init(double alt, double emin, double emax) {
// initialize the model parameters assuming altitude "alt" (km)

    m_version = 18052005; // change this version number if you change any of the parameters below!

    if(alt < 0.){
        alt = 426.79;
    }
    else if(alt < 300.){
        std::cout << "WARNING: Earth model invalid below 300 km altitude." << std::endl;
    }
    if(emin < 10.){
        std::cout << "WARNING: Earth model invalid below 10 MeV." << std::endl;
    }
    if(emax > 350000.){
        std::cout << "WARNING: Earth model only an extrapolation above 350 GeV." << std::endl;
    }
    if(emin > emax){
        std::cout << "ERROR in Earth model input parameters: Emax must be larger than Emin." 
		  << std::endl;
        throw std::invalid_argument("Earth: ERROR in Earth model input parameters: Emax must be larger than Emin.");
    }

// Earth radius (km)
    double re = 6378.14;
// Orbit radius (km)
    double ro = re + alt;


// Geometrical Horizon ZA (deg)
    m_a[1] = 90. + 180./3.1415926 * acos(re/ro);
// altitude scale factor
    m_a[2] = 0.; // set futher below, when m_a[18] is defined
// thickness of atmosphere (km)
    m_a[3] = 100.; 
// distance spacecraft - horizon (km)
    m_a[4] = ro * sqrt(1. - re*re/(ro*ro));
// normalization of g_0(E) (cm^(-2)s^(-1)sr^(-1)MeV^(-1))
    m_a[5] = exp(-0.4716); 
// index of g_0(E)  
    m_a[6] = -2.088;
// cutoff of g_0(E) (MeV) 
    m_a[7] = 3E4; // (30 GeV) 
// normalization of g_1(E) (cm^(-2)s^(-1)sr^(-1)MeV^(-1))
    m_a[8] = exp(-1.036); 
// index of g_1(E) 
    m_a[9] = -1.811; 
// cutoff of g_1(E) (MeV)
    m_a[10] = 2914.; 
// normalization of \sigma_(AZ)(\theta) (deg) 
    m_a[11] = 76.9; 
// mean of \sigma_(AZ)(\theta) (deg)
    m_a[12] = 98.6;
// std. dev. of \sigma_(AZ)(\theta) (deg)
    m_a[13] = 13.8;
// peak position of azimuthal profile, 180^\circ = West (deg)
    m_a[14] = 180.;
// normalization of central flux  (cm^(-2)s^(-1)sr^(-1)MeV^(-1))
    m_a[15] = exp(-0.06731);
// index of central flux 
    m_a[16] = -2.512; 
// cutoff of central flux (MeV)
    m_a[17] = 3000.; 
// reference geometrical horizon position (deg)
    m_a[18] = 110.4;
    m_a[2] = (180. - m_a[1])/(180. - m_a[18]);
// altitude (km)
    m_a[19] = alt;
// reference altitude (km)
    m_a[20] = 426.79;
// minimum energy (MeV)
    m_a[21] = emin;
// maximum energy (MeV)
    m_a[22] = emax;
// integral flux over whole sky (cm^(-2)s^(-1))
    m_ftot = ff(200);

//     for(i=1;i<=22;i++){
//       std::cout << i << ": " << m_a[i] << " ";
//     }
//     std::cout << std::endl;

    return;
}
       
double Earth::pp(double t, double p, double e) const {
// the probability density normalized to give 1 when integrated over the
//  whole sky
    return f(t,p,e) / m_ftot;
}

double Earth::qn(double e) const {
// the normalization of the envelope function
    return pp(tp(),180.,e)/pow(e,-1.5);
}

double Earth::q(double e) const {
// the envelope function
    return qn(m_a[21])*pow(e,-1.5);
}

double Earth::qq(double e) const {
// the integral of the envelope function
    return 4.* 3.141592 * qn(m_a[21]) * 2. * (pow(m_a[21],-0.5) - pow(e,-0.5));
}

double Earth::qqinv(double x) const {
// the inverse of the integral of q
    return pow((pow(m_a[21],-0.5) - x/(8.*3.141592*qn(m_a[21]))),-2.);
}


void Earth::earth(double &t, double &p, double &e) const {
// return zenith angle (deg), azimuth (deg) and energy (MeV)
    double e0,t0,p0,aqq;
    double dummy;
    double radpdg;
    int n,count;
    n = 1;  
    radpdg = 3.1415926/180.;
    aqq = qq(m_a[22]);
    count = 0;
    dummy = 0.;
    while(count < n){
        dummy = dummy + 1.;
        double r1=CLHEP::RandFlat::shoot(),
               r2=CLHEP::RandFlat::shoot(),
               r3=CLHEP::RandFlat::shoot(),
               r4=CLHEP::RandFlat::shoot();
        
        e0 = qqinv(r1*aqq);
// uniform 2D coordinates in one hemisphere,
// we know that pp is zero at t0 < 90.
        t0 = 180. - acos(r2)/radpdg;
        p0 = r3 * 360.;
        if(r4*q(e0) < pp(t0,p0,e0)){
            t = t0;
            p = p0;
            e = e0;
            ++count;
        }
    }
    return;
}



ISpectrumFactory &Magrathea() { // a.k.a. EarthFactory, see http://www.bbc.co.uk/dna/h2g2/A105265
   static SpectrumFactory<Earth> myFactory;
   return myFactory;
}

Earth::Earth(const std::string &paramString) {

   std::vector<std::string> params;
   facilities::Util::stringTokenize(paramString, ", ", params);

   m_alt = std::atof(params[0].c_str());
   m_emin = std::atof(params[1].c_str());
   m_emax = std::atof(params[2].c_str());

   init(m_alt, m_emin, m_emax);

   std::cerr << "Earth created. Total flux = " 
	     << m_ftot << " cm^-2 s^-1 " << " between " 
	     << m_emin << " MeV and "
	     << m_emax << " MeV." << std::endl;

   m_eCalled = false;

}

double Earth::flux(double time) const { // argument is the mission elapsed time (s)
    return m_ftot * 1E4 /solidAngle(); // (m^-2 / s/ sr)
}

double Earth::solidAngle() const {
  double rval;
  rval = 2.*M_PI*(1. - cos( M_PI/180. * ( 180.-tp()+5*sigma() ) ) );
  return rval;
}


double Earth::energy(double time) {
    (void)(time); // the Earth is not variable in this version
    earth(m_t, m_p, m_e);
    m_eCalled = true;
    return (m_e);
}

std::pair<double, double> Earth::dir(double energy) {
    double costheta, phi;
    if(energy != m_e || !m_eCalled){
	std::cerr << "ERROR in routine calling Earth: need to call energy() before dir()." << std::endl;
        throw std::runtime_error("Earth: ERROR in routine calling Earth: need to call energy() before dir().");
    }
    m_eCalled = false;
// Assume ZA AZ coordinates by default.
    costheta = cos(m_t/180.*M_PI);
    phi = m_p/180.*M_PI; // AZ is 0 in the East and 90. deg in the North when looking at the Earth
    return std::make_pair(costheta, phi);
}

std::pair<double, std::pair<double, double> > Earth::photon(){
    earth(m_t, m_p, m_e);
    std::pair<double, double> theDir = std::make_pair(cos(m_t/180.*M_PI), m_p/180.*M_PI);
    return std::make_pair(m_e, theDir);
}
