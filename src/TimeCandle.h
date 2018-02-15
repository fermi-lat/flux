/**
 * @file TimeCandle.h
 * @brief Declaration of class TImeCandle.cxx: a source that ticks 

 * $Header: /nfs/slac/g/glast/ground/cvs/flux/src/TimeCandle.h,v 1.5 2006/11/15 20:57:42 burnett Exp $
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
* $Header: /nfs/slac/g/glast/ground/cvs/flux/src/TimeCandle.h,v 1.5 2006/11/15 20:57:42 burnett Exp $
*/


class TimeCandle : public Spectrum {
public: 
    TimeCandle(const std::string& params);
    
    TimeCandle();
    virtual const char* particleName()const;
    virtual std::string title()const;
        
    virtual std::pair<double,double> dir(double){
        return std::make_pair(1.0,0.0);
    }     
    
    double energy( double time);
    
    double interval (double time);
private:
    float parseParamList(std::string input, int index, float default_value);
    double m_period; ///< period 
    double m_offset; ///< offset in seconds after second boundary to start (-1: ignore)
    std::string m_name;	// particle name to generate ("P", "gamma", ...)
    bool m_first;
};

#endif
