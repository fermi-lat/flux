// $Header: /nfs/slac/g/glast/ground/cvs/flux/src/TimeCandle.h,v 1.1.1.1 2003/07/29 18:22:19 burnett Exp $
#ifndef TimeCandle_H
#define TimeCandle_H

#include "flux/Spectrum.h"
#include <string>

/** 
* \class TimeCandle
*
* @brief Spectrum: base class for energy spectrum objects
* TimeCandle: define a particle with a constant time of arrival.

  a convenient Spectrum : a single particle at a constant incremental time, 
* 
* $Header: /nfs/slac/g/glast/ground/cvs/flux/src/TimeCandle.h,v 1.1.1.1 2003/07/29 18:22:19 burnett Exp $
*/


class TimeCandle : public Spectrum {
public: 
    TimeCandle(const std::string& params);
    
    TimeCandle();
    //void setPosition ( float /*lat*/, float /*lon*/ ){}
    //virtual double calculate_rate(double old_rate);
    virtual float  operator()(float f);
    virtual const char* particleName()const;
    virtual std::string title()const;
        
    virtual std::pair<float,float> dir(float)const{
        return std::make_pair<float,float>(1.0,0.0);
    }     
    
    virtual std::pair<double,double> dir(double energy){return std::make_pair<double,double>(1.0,0.0);}
    
    double energy( double time){  return 0.;}
    
    double interval (double time)
    {        
        return m_T0;
    }
private:
    float parseParamList(std::string input, int index);
    double m_T0; //how many seconds pass before each next incoming particle
    std::string m_name;	// particle name to generate ("P", "gamma", ...)
};

#endif
