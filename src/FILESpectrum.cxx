/** 
* @file FILESpectrum.cxx
* @brief Implementation of FILESpectrum
*
*  $Header: /nfs/slac/g/glast/ground/cvs/flux/src/FILESpectrum.cxx,v 1.4 2004/12/21 03:46:40 burnett Exp $
*/

#include "FILESpectrum.h"

#include <cmath>
#include <algorithm>
#include <functional>
#include <fstream>
#include <iostream>

#include "flux/SpectrumFactory.h"
#include "facilities/Util.h"

static SpectrumFactory<FILESpectrum> factory;
const ISpectrumFactory& FILESpectrumFactory = factory;

FILESpectrum::FILESpectrum(const std::string& params) : m_inLog(false)
  //  :m_inGeV(true),m_inLog(false),m_particle_name("p")
{
  std::vector<efpair> temp_vector = initialize(params);
        
  double total_flux = 0.0;
  std::pair<double,double> ef;
  std::vector<std::pair<double,double> >::iterator it; 
  for(it = temp_vector.begin(); it != temp_vector.end(); it++)
    {
      ef = (*it);
      double factor;
      
      if(it == temp_vector.begin() && (it+1) != temp_vector.end() )
	factor = ( (it+1)->first - ef.first );
      else if( it == temp_vector.begin() && (it+1) == temp_vector.end() )
	factor = 1;
      else if( (it+1) != temp_vector.end() )
	factor = (((it+1)->first - (it-1)->first)/2);
      else if( (it+1) == temp_vector.end() )
	factor = ( ef.first - (it-1)->first );
      
//      std::cout<<factor<<" "<<ef.second<<std::endl;
      ef.second *= factor;
      total_flux += ef.second;
      
      ef.second = total_flux;
      
      integ_flux.push_back(ef);
    }
  m_fileflux = total_flux;
}


double FILESpectrum::flux() const
{
    /// Purpose: Calculate flux for the current position
    if(m_flux==0) 
     return m_fileflux;
    else
     return m_flux;
}

double FILESpectrum::flux (double time ) const{
    return flux();
}


float FILESpectrum::operator() (float r)
{
    /// Purpose: sample a single particle energy from the spectrum
    double target_flux = r * m_fileflux;
    
    std::vector<efpair>::const_iterator i;
    
    std::pair<double,double> previous;
    
    i = integ_flux.begin();
    previous = (*i);
    
    for(i = integ_flux.begin(); i != integ_flux.end(); i++)
    {
        if((*i).second >= target_flux)
            break;
        previous = (*i);
    }
    
    // Use linear interpolation between bins
    double m = ( (*i).first - previous.first ) / ( (*i).second - previous.second );
    double b = (*i).first - m * (*i).second;
    float scale = 1.;
    double raw_e = m * target_flux + b;
//    if(!m_inGeV) scale = 0.001;
    if(m_inLog)
      {
        return scale*pow(10., raw_e);
      }
    else
      {
        return scale*raw_e;
      }
}


std::string FILESpectrum::title() const
{
    return "FILESpectrum";
}

const char * FILESpectrum::particleName() const
{
    return m_particle_name.c_str();
}

std::vector<std::pair<double,double> > FILESpectrum::initialize(const std::string& params)
{
  //initialize some variables:
  const char* flux_root = ::getenv("FLUXROOT");
  std::string doc_path= (flux_root? std::string(flux_root)+"/" : "");
  std::string fileName = doc_path;
  std::ifstream input_file;
  
  std::vector<std::string> output;
  facilities::Util::stringTokenize(params,",",output);
  
  for(std::vector<std::string>::const_iterator it = output.begin();it!=output.end();it++)
    {
      if(! input_file.is_open())
	{
	  fileName.append((*it));
	  input_file.open(fileName.c_str());
	} else if(*it == "log") 
	    {
	      m_inLog = true; 
	    } 
    }
      
  if(! input_file.is_open())
    {
      std::cerr << "ERROR:  Unable to open:  " << fileName.c_str() << std::endl;
      throw(std::string("ERROR:  Unable to open: "+ fileName ));
    }
  else 
    {
      return readFile(input_file);
    }

}

std::vector<std::pair<double,double> > FILESpectrum::readFile(std::ifstream& input_file)
{
  char buffer[256];
  std::vector<std::pair<double,double> > temp_vector;

  while(input_file.getline(buffer,256,'\n'))
    {
      std::string line(buffer);
      std::vector<std::string> entries;
      if(line.find('%')!= std::string::npos)
      	{
	  //This is taken to be a comment line
	  continue;
      	}
      facilities::Util::stringTokenize(line," \t",entries);
      if(entries.back() == "") entries.pop_back(); //this should be removed after a fix in stringTokenize
      int size = entries.size();
      if(size == 2)
	{
	  temp_vector.push_back(std::make_pair<double,double>(atof(entries[0].c_str()),atof(entries[1].c_str())));
	} 
      else if(size == 3)
	{
	  std::cerr<<"Line contained more than 2 columns"<<std::endl;
        }
      
    }
  return temp_vector;
}
