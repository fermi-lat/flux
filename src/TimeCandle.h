/**
 * @file TimeCandle.h
 * @brief Declaration of class TImeCandle.cxx: a source that ticks 

 * $Header$
 */

#ifndef flux_TimeCandle_H
#define flux_TimeCandle_H

#include "flux/Spectrum.h"
#include <string>

/** 
* \class TimeCandle
*
* @brief Spectrum: base class for energy spectrum objects
* TimeCandle: define a particle with a constant time of arrival.

  a convenient Spectrum : a single particle at a constant incremental time, 
  @author: S. Robinson
* 
* $Header: /nfs/slac/g/glast/ground/cvs/flux/src/TimeCandle.h,v 1.3 2005/03/18 04:38:17 burnett Exp $
*/


class TimeCandle : public Spectrum {
public: 
    TimeCandle(const std::string& params);
    
    TimeCandle();
    virtual const char* particleName()const;
    virtual std::string title()const;
        
    virtual std::pair<double,double> dir(double){
        return std::make_pair<float,float>(1.0,0.0);
    }     
    
    double energy( double time);
    
    double interval (double time);
private:
    float parseParamList(std::string input, int index);
    double m_T0; //how many seconds pass before each next incoming particle
    std::string m_name;	// particle name to generate ("P", "gamma", ...)
    bool m_first;
};

#endif
