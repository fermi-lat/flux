// Spectrum.cxx: implementation of the Spectrum class.
//
//////////////////////////////////////////////////////////////////////

#include "flux/AGNSpectrum.h"
#include <cmath>

// CLHEP
#include "CLHEP/Random/RandomEngine.h"
#include "CLHEP/Random/RandGeneral.h"
#include "CLHEP/Random/RandExponential.h"
#include "CLHEP/Random/RandFlat.h"

// this is needed to include in the executable or dll
#include "flux/SpectrumFactory.h"

static SpectrumFactory<AGNSpectrum> factory;
const ISpectrumFactory& AGNSpectrumFactory = factory;

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

AGNSpectrum::AGNSpectrum(const std::string& params){
    std::vector<float> input1;
    std::vector<float>& inputs = input1;
    parseParamList(params, inputs);
    init(inputs);
    setInitialFlareStates();
}

void AGNSpectrum::init(std::vector<float> params){

    m_l =  params.size()>0? params[0]: 0.0f;
    m_b =  params.size()>1? params[1]: 0.0f;
    m_flux =  params.size()>2? params[2]: 10.0f;
    m_index =  params.size()>3? params[3]: 2.0f;
    m_flareMult =  params.size()>4? params[4]: 5.0f;
    m_flareAdd =  params.size()>5? params[5]: 0.01f;
    m_flareDuty =  params.size()>6? params[6]: 0.03f;
    m_flarePeriod =  params.size()>7? params[7]: 1000.0f;
    m_Emin =  params.size()>8? params[8]: 0.1f;
    m_Emax =  params.size()>9? params[9]: 100.0f;
    

}

AGNSpectrum::AGNSpectrum(float l,float b,float flux,float index,float flareMult,float flareAdd,float flareDuty,float flarePeriod,float Emin,float Emax):
m_l(l),m_b(b),m_flux(flux),m_index(index),m_flareMult(flareMult),m_flareAdd(flareAdd),m_flareDuty(flareDuty),m_flarePeriod(flarePeriod),m_Emin(Emin),m_Emax(Emax)
{setInitialFlareStates();}

AGNSpectrum::~AGNSpectrum(){}

float AGNSpectrum::operator()(float r)const{
    double index;
    //double E0=0.1;
    //double emax=100.;
    double rand=RandFlat::shoot(1.);

    if(m_flaring){
        index=m_index+m_flareAdd;
    }else{
        index=m_index;
    }
    float x = 1 - exp((1-index)*log(m_Emax/0.5/*m_Emin*/));
    return 0.5/*m_Emin*/*exp(log(1-x*rand)/(1-index));
}

double AGNSpectrum::flux (double time ) const {
    if(!m_flaring)return m_flux;
    return m_flux*m_flareMult;
}

double AGNSpectrum::solidAngle( )const
{
    return 1.0; // flag that doesn't calculate
}


double AGNSpectrum::energy(double time)
{
    // default implementation, which works for other Spectrum objects
    return (*this)(time);
}


std::pair<double,double> AGNSpectrum::dir(double energy)
{
    return std::make_pair<double,double>(m_l,m_b);
}

void AGNSpectrum::parseParamList(std::string input, std::vector<float>& output) const
{

    int i=0;
    for(;!input.empty() && i!=std::string::npos;){
        float f = ::atof( input.c_str() );
        output.push_back(f);
        i=input.find_first_of(",");
        input= input.substr(i+1);
    } 
}


double AGNSpectrum::interval (double time)
{
    double r;
    double calculatedInterval=1E8;

	double ECutoffMult=pow(1.0/m_Emin,-1.0*(m_index-1));
    r = (solidAngle()*flux(time)*ECutoffMult*60000./*EventSource::totalArea()*/);

    if (r == 0){ return -1.;
    }else{  
        double p = RandFlat::shoot(1.);
        //std::cout << "interval from the AGN:= " << (-1.)*(log(1.-p))/r << std::endl;
        calculatedInterval = (-1.)*(log(1.-p))/r;
    }
    //ok, this rate only applies until the next state flip, then the interval must be calculated again:
    if(time + calculatedInterval > m_nextFlipTime){
        flipState();
        reCalcNextFlip();
        return interval(time + calculatedInterval);
    }
    return calculatedInterval;
}

void AGNSpectrum::setInitialFlareStates(){
    m_flaring=false;
    m_nextFlipTime=0.;
    reCalcNextFlip();
}

void AGNSpectrum::reCalcNextFlip(){
    double howLong; //how long until the next time flip (flaring/quiescent).
    if(m_flaring){
        howLong=m_flarePeriod;
    }else{
        howLong=m_flarePeriod*((1.-m_flareDuty)/(m_flareDuty));

    }
    m_nextFlipTime+=howLong;
}
