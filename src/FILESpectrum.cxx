/** 
* @file FILESpectrum.cxx
* @brief Implementation of FILESpectrum
*
*  $Header: /nfs/slac/g/glast/ground/cvs/FluxSvc/src/FILESpectrum.cxx,v 1.8 2002/11/20 20:09:54 srobinsn Exp $
*/

#include "FILESpectrum.h"

#include <cmath>
#include <algorithm>
#include <functional>
#include <fstream>
#include <iostream>

#include "flux/SpectrumFactory.h"

static SpectrumFactory<FILESpectrum> factory;
const ISpectrumFactory& FILESpectrumFactory = factory;

FILESpectrum::FILESpectrum(const std::string& params)
{
    m_particle_name = "p";
    // TODO: have params be an index 
    
    if(params.empty())
        initialization_document = "sources/glast_smin_flx.txt";
    else
        initialization_document = params.c_str();
    
    // construct filename, assuming data files in folder <packageroot>/CREME
    std::string fileName = "";
    const char* flux_root = ::getenv("FLUXROOT");
    std::string doc_path= (flux_root? std::string(flux_root)+"/" : "");
    fileName = doc_path+initialization_document;
    std::ifstream input_file;
    input_file.open(fileName.c_str());
    
    if(! input_file.is_open())
    {
        std::cerr << "ERROR:  Unable to open:  " << fileName.c_str() << std::endl;
        throw(std::string("ERROR:  Unable to open: "+ fileName ));
    }
    else
    {
        const int buffer_size = 256;
        char buffer[buffer_size];
        input_file.getline(buffer,buffer_size,'\n');
        std::pair<double,double> ef;
        std::vector<efpair> temp_vector;
        double total_flux = 0;
        {/* Skip whitespace */
            char temp = input_file.peek();
            while(temp == ' ' || temp == '\n' || temp == '\t')
            {
                temp = input_file.get();
                temp = input_file.peek();
            }
        }
        while('.' == input_file.peek() || isdigit(input_file.peek()))
        {
            input_file >> ef.first;
            input_file >> ef.second;
            
            {/* Convert units */
                ef.first /= 1000;  // MeV -> GeV
                ef.second *= 1000; // particles/m2-s-sr-MeV -> particles/m2-s-sr-GeV
            }
            
            temp_vector.push_back(ef);
            
            {/* Skip whitespace */
                char temp = input_file.peek();
                while(temp == ' ' || temp == '\n' || temp == '\t')
                {
                    temp = input_file.get();
                    temp = input_file.peek();
                }
            }
        }
        
        std::vector<efpair>::iterator it; 
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
            
            ef.second *= factor;
            total_flux += ef.second;
            
            ef.second = total_flux;
            
            integ_flux.push_back(ef);
        }
        m_flux = total_flux;
    }
}


double FILESpectrum::flux() const
{
    /// Purpose: Calculate flux for the current position
    return m_flux;
}

double FILESpectrum::flux (double time ) const{
    return flux();
}


float FILESpectrum::operator() (float r)const
{
    /// Purpose: sample a single particle energy from the spectrum
    double target_flux = r * m_flux;
    
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
    
    return m * target_flux + b;
}


std::string FILESpectrum::title() const
{
    return "FILESpectrum";
}

const char * FILESpectrum::particleName() const
{
    return m_particle_name.c_str();
}

