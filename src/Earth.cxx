/**
 * @file Earth.cxx
 * @brief A phenomenological model of the Earth based on EGRET measurements
 * @author D. Petry
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/terrestrialSources/Earth/src/Earth.cxx,v 1.3 2004/12/14 23:37:57 petry Exp $
 */

#include <iostream>

#include <cmath>
#include <cstdlib>

#include "facilities/Util.h"

#include "flux/SpectrumFactory.h"
#include "flux/EventSource.h"

#include "Earth.h"

static SpectrumFactory<Earth> factory;
const ISpectrumFactory& EarthFactory = factory;


double Earth::tp(double ee) const {
// Peak position
    double e;
    if(ee < 35.){
        e = 35.;
    } else{
        e = ee;
    }
    return m_a[1] + m_a[2]* pow(e,m_a[3]) * m_a[21]/m_a[1] ;
}

double Earth::sigma(double ee) const {
// Peak Width corrected for EGRET PSF, orbital decay, orbital interpolation
    double e;
    if(ee < 35.){
        e = 35.;
    } else{
        e = ee;
    }
    return m_a[4]*pow(e,m_a[5]) * m_a[21]/m_a[1];
}

double Earth::sigmaUncorrected(double ee) const {
// Peak Width as measured by EGRET
    double e;
    if(ee < 35.){
        e = 35.;
    } else{
        e = ee;
    }
    return m_a[26]*pow(e,m_a[27]) * m_a[21]/m_a[1];
}

double Earth::t0(double ee) const {
// Transition point from Gaussian to Exponential 
    double e;
    if(ee < 35.){
        e = 35.;
    } 
    else{
        e = ee;   
    }
    if(e > 60.){
        return tp(e) + ( m_a[6] + sqrt(log10(e-m_a[7])) ) * sigma(e) * m_a[21]/m_a[1];
    } else{
        return tp(e) * m_a[21]/m_a[1];
    }
}

double Earth::g0(double e) const {
// Spectrum of the constant component independent of AZ
    return m_a[8]*pow(e,m_a[9])*exp(-e/m_a[10]);
}

double Earth::g1(double e) const {
// Spectrum of the component dependent on AZ
    return m_a[11]*pow(e,m_a[12])*exp(-e/m_a[13]);
}

double Earth::sigmaz(double t) const {
// Std Dev of the AZ dependent component
    return m_a[14] * exp(-0.5*pow(((t-m_a[15])/m_a[16]),2.));
}

double Earth::f0(double t, double p, double e) const {
// the part zt ZA < t0
    return g0(e) + g1(e) * exp(-0.5*pow(((p-m_a[17])/sigmaz(t)),2.));
}

double Earth::normcorr(double e) const {
// correction of the normalization of fa(t,p,e)
// necessary because we are replacing the measured peak width 
// by sigma(e) which is corrected for the EGRET PSF and effects
// from orbital decay and orbital interpolation.
  return sigmaUncorrected(e)/sigma(e);
}

double Earth::fa(double t, double p, double e) const {
    if(t <= t0(e)  &&  t > t0(e)-5.*sigma(e)){
        return normcorr(e)*exp(-0.5*pow(((t-tp(e))/sigma(e)),2.))*f0(t,p,e);
    } else{
        return 0.;
    }
}

double Earth::b(double e) const {
// central flux
    return m_a[18]*pow(e,m_a[19])*exp(-e/m_a[20]);
}

double Earth::c(double e, double p) const {
    return (log(fa(t0(e),p,e))-log(b(e)))/(t0(e)-180.);
}

double Earth::c0(double e, double p) const {
    return log(fa(t0(e),p,e))-c(e,p)*t0(e);
}

double Earth::fb(double t, double p, double e) const {
    if(t > t0(e)){
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
// Integral of earth flux over all sky from emin to emax
//  using isteps integration steps
    int i, j, k;
    double theint, t, p, dt, dp, dloge, domega, dtheint;
    double radpdg;
    double steps, logemin, logemax;
    double f1, f2, e1, e2, t1, t2;
    double norm, alpha;
// check if precomputed value available
    theint = 0.;
    if(m_a[24] == 10. && m_a[25] == 10000.){
        if(isteps == 100){
            if(m_a[22] == 565.){
		if(m_version == 14122004){
		    std::cout << "Earth: Using precomputed value for integral flux." << std::endl;
		    theint = 0.181899;
		}
            }
        }
    }
    if(theint == 0.){
        radpdg = M_PI/180.;  
        steps = isteps/1.;
        logemin = log(m_a[24]);
        logemax = log(m_a[25]);
        dp = 360./steps;
        dloge =(logemax-logemin)/steps;
        dt = 180./steps;

        for(i=0;i<isteps;i++){
            
//            if(i/100. == i/100){
//                std::cout << i <<" of " << isteps << std::endl;
//            }

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
            norm = f1/pow(e1,-alpha);
            dtheint = norm * ( 
                pow(e2,-alpha+1.)/(-alpha+1.)
                -  pow(e1,-alpha+1.)/(-alpha+1.)
                );
            theint = theint + dtheint;
        }
    }
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
    double norm, alpha;
    double ringFraction;

    theint = 0.;
    omega = 0.;
    if(theint == 0.){
        radpdg = M_PI/180.;  
        steps = isteps/1.;
        logemin = log(m_a[24]);
        logemax = log(m_a[25]);
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
            norm = f1/pow(e1,-alpha);
            dtheint = norm * ( 
                pow(e2,-alpha+1.)/(-alpha+1.)
                -  pow(e1,-alpha+1.)/(-alpha+1.)
                );
            theint = theint + dtheint;
        }
    }
    omega = omega/steps;
    std::cout << "omega = " << omega << std::endl;
    return theint;
}


void Earth::init(double alt, double emin, double emax){
// initialize the model parameters assuming altitude "alt" (km)

    m_version = 14122004; // change this version number if you change any of the parameters below!

    double re,ro;
// Earth radius (km)
    re = 6378.14;
// Orbit radius (km)
    ro = re + alt;

    if(alt < 0.){
        alt = 426.79;
    }
    else if(alt < 300.){
	std::cout << "WARNING: Altitude = " << alt 
		  << ". Earth model invalid below 300 km altitude." << std::endl;
    }
    if(emin < 10.){
        std::cout << "WARNING: Emin = " << emin
		  << ". Earth model invalid below 10 MeV." << std::endl;
    }
    if(emax > 10000.){
        std::cout << "WARNING: Emax = " << emax
		  << ". Earth model invalid above 10 GeV." << std::endl;
    }
    if(emin > emax){
        std::cout << "ERROR in Earth model input parameters: Emin = " << emin
		  << ", Emax = " << emax << ". Emin must be larger than Emax." << std::endl;
	exit(1);
    }

// Geometrical Horizon ZA (deg)
    m_a[1] = 90. + 180./3.1415926 * acos(re/ro);
// normalization of \theta_(peak)(E) (deg)
    m_a[2] = exp(3.6835); 
// index of \theta_(peak)(E)
    m_a[3] = -0.4830; 
// normalization of \sigma(E) (deg)
    m_a[4] = exp(3.011); 
// index of \sigma(E)
    m_a[5] = -0.350; 
// energy dependence of \theta_(0)
    m_a[6] = -0.569;  
// energy dependence of \theta_(0) (MeV)
    m_a[7] = 57.45; 
// normalization of g_0(E) (cm^(-2)s^(-1)sr^(-1)MeV^(-1))
    m_a[8] = exp(-0.8219); 
// index of g_0(E)  
    m_a[9] = -2.000;
// cutoff of g_0(E) (MeV) 
    m_a[10] = 2514.; 
// normalization of g_1(E) (cm^(-2)s^(-1)sr^(-1)MeV^(-1))
    m_a[11] = exp(-1.036); 
// index of g_1(E) 
    m_a[12] = -1.811; 
// cutoff of g_1(E) (MeV)
    m_a[13] = 2914.; 
// normalization of \sigma_(AZ)(\theta) (deg) 
    m_a[14] = 76.9; 
// mean of \sigma_(AZ)(\theta) (deg)
    m_a[15] = 98.6;
// std. dev. of \sigma_(AZ)(\theta) (deg)
    m_a[16] = 13.8;
// peak position of azimuthal profile, 180^\circ = West (deg)
    m_a[17] = 180.;
// normalization of central flux  (cm^(-2)s^(-1)sr^(-1)MeV^(-1))
    m_a[18] = exp(-0.06731);
// index of central flux 
    m_a[19] = -2.512; 
// cutoff of central flux (MeV)
    m_a[20] = 3000.; 
// reference geometrical horizon position
    m_a[21] = 110.4;
// altitude (km)
    m_a[22] = alt;
// reference altitude (km)
    m_a[23] = 426.79;
// minimum energy (MeV)
    m_a[24] = emin;
// maximum energy (MeV)
    m_a[25] = emax;
// normalization of uncorrected \sigma(E) (deg)
    m_a[26] = exp(3.0385); 
// index of uncorrected \sigma(E)
    m_a[27] = -0.2825; 
// integral flux over whole sky (cm^(-2)s^(-1))
    m_ftot = ff(100);

//    for(i=1;i<=27;i++){
//        std::cout << m_a[i] << " ";
//    }
//    std::cout << std::endl;
    return;
}
      
double Earth::pp(double t, double p, double e) const {
// the probability density normalized to give 1 when integrated over the
//  whole sky
    return f(t,p,e) / m_ftot;
}

double Earth::qn(double e) const {
// the normalization of the envelope function
    return pp(tp(e),180.,e)/pow(e,-1.5);
}

double Earth::q(double e) const {
// the envelope function
    return qn(m_a[24])*pow(e,-1.5);
}

double Earth::qq(double e) const {
// the integral of the envelope function
    return 4.* 3.141592 * qn(m_a[24]) * 2. * (pow(m_a[24],-0.5) - pow(e,-0.5));
}

double Earth::qqinv(double x) const {
// the inverse of the integral of q
    return pow((pow(m_a[24],-0.5) - x/(8.*3.141592*qn(m_a[24]))),-2.);
}


void Earth::earth(double &t, double &p, double &e) const {
// return zenith angle (deg), azimuth (deg) and energy (MeV)
    double r1,r2,r3,r4;
    double e0,t0,p0,aqq;
    double dummy;
    double radpdg;
    int n,count;
    n = 1;  
    radpdg = 3.1415926/180.;
    aqq = qq(m_a[25]);
    count = 0;
    dummy = 0.;
    while(count < n){
        dummy = dummy + 1.;
	
        r1 = rand()/(RAND_MAX+1.0);
        r2 = rand()/(RAND_MAX+1.0);
        r3 = rand()/(RAND_MAX+1.0);
        r4 = rand()/(RAND_MAX+1.0);
        
        e0 = qqinv(r1*aqq);
// uniform 2D coordinates in one hemisphere,
// we know that pp is zero at t0 < 90.
        t0 = 180. - acos(r2)/radpdg;
        p0 = r3 * 360.;
        if(r4*q(e0) < pp(t0,p0,e0)){
            t = t0;
            p = p0;
            e = e0;
            count = count + 1;
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
    (void)(time); // no variability in this version
    return m_ftot * 1E4; // (m^-2 s^-1)
}

double Earth::solidAngle() const {
  double rval;
  double e = 35.; // we return the largest possible solid angle, i.e. the one at 35 MeV
  rval = 2.*M_PI*(1. - cos( M_PI/180. * ( 180.-tp(e)-5*sigma(e) ) ) );
  return rval;
}

double Earth::interval(double time) { // argument is the mission elapsed time (s)
   double theInterval;
   double rate = flux(time)*EventSource::totalArea();
   double xi = rand()/(RAND_MAX+1.0);
   if(rate != 0.){
       theInterval = -log(1. - xi)/rate;
   }
   else {
       theInterval = 5E17;
   }
   return theInterval;
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
	exit(1);
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
