/** 
* @file GalPulsars.cxx
* @brief definition of GalPulsars
*
*/

#include "flux/Spectrum.h"
#include "flux/SpectrumFactory.h"
#include "flux/EventSource.h"
#include "astro/EarthOrbit.h"
#include "CLHEP/Random/RandFlat.h"
#include "facilities/Util.h"
#include <string>
#include <utility>
#include <algorithm>
#include <vector>
#include <fstream>

#define GAL_SWITCH


class GalPulsars : public Spectrum
{
public:

   GalPulsars(const std::string& params);

   std::pair<double,double> dir(double energy);
   float operator()(float xi) const;

   double interval(double time);
   double energy(double time);

   virtual std::string title() const{return "GalPulsars";}
   virtual const char * particleName() const {return "gamma";}
   inline  const char * nameOf() const {return "GalPulsars";}

private:

   void updateIntervals(double current_time, double time_decrement);
   double period(double time, int pulsarIndex) const;
   double power_law( double r, double e1, double e2, double gamma) const;
   void initLightCurve(void);

   std::vector< std::vector<double> > m_lightCurve;
   std::vector< std::string > m_name;
   std::vector<double> m_spectralIndex;
   std::vector<double> m_highCutoff;
   std::vector<double> m_lowCutoff;
   std::vector<double> m_flux;
   std::vector<double> m_freq_dot;
   std::vector<double> m_freq;
   std::vector<double> m_lat;
   std::vector<double> m_lon;
   std::vector<double> m_t0;

   std::vector<double> m_interval;
   int m_changed;
};

static SpectrumFactory<GalPulsars> factory;
const ISpectrumFactory& GalPulsarsFactory = factory;

GalPulsars::GalPulsars(const std::string& paramString)
{
   std::vector<std::string> params;
   facilities::Util::stringTokenize(paramString, ", ", params);

   std::ifstream input_file;
   input_file.open(params[0].c_str(), std::ios::in);

   if(!input_file.is_open())
   {
      std::cerr << "Error:  Unable to read input file" << std::endl;
      return;
   }

   // Skip first line containing header information
   char buffer[1024];
   input_file.getline(buffer,1024);

   while(!input_file.eof())
   {
      input_file.getline(buffer,1024,'\t');
      if(buffer[0] != 'G')
         break;

      m_name.push_back(buffer);



      input_file.getline(buffer,1024,'\t'); 
      m_lat.push_back(std::atof(buffer));

      input_file.getline(buffer,1024,'\t'); 
      m_lon.push_back(std::atof(buffer));

      input_file.getline(buffer,1024,'\t'); 
      m_freq.push_back(1./std::atof(buffer));

      input_file.getline(buffer,1024,'\t'); 
      m_freq_dot.push_back(-std::atof(buffer)*m_freq.back()*m_freq.back());

      input_file.getline(buffer,1024,'\t');
      m_flux.push_back(std::atof(buffer));

      input_file.getline(buffer,1024,'\t');
      m_highCutoff.push_back(std::atof(buffer));

      m_lowCutoff.push_back(0.03);

      input_file.getline(buffer,1024,'\t');
      m_spectralIndex.push_back(std::atof(buffer));

      m_t0.push_back(0.0);

      std::vector<double> temp_lc;
      for(int i = 0; i < 20; i++)
      {
         input_file.getline(buffer,1024,'\t');
         temp_lc.push_back(std::atof(buffer));
      }
      input_file.getline(buffer,1024,'\n');
      temp_lc.push_back(std::atof(buffer));

      m_lightCurve.push_back(temp_lc);
      temp_lc.clear();
   }
   input_file.close();

   m_changed = -1;
   m_interval.resize(m_flux.size());

   initLightCurve();
}


double GalPulsars::interval(double current_time) 
{
   if(m_changed == -1)
      updateIntervals(current_time,0.);

   // Find the index of the shortest interval
   m_changed = 0;
   for(unsigned int i = 1; i < m_interval.size(); i++)
   {
      if(m_interval[i] < m_interval[m_changed]) 
         m_changed = i;
   }

   // Store the shortest interval before calling the update
   double shortestinterval = m_interval[m_changed];
   
   // Update the intervals
   updateIntervals(current_time, shortestinterval);

   return shortestinterval;
}

void GalPulsars::updateIntervals(double current_time, double time_decrement)
{
   for(unsigned int i = 0; i < m_interval.size(); i++)
   {
      if(m_changed == i || m_changed == -1)
      {
         astro::EarthOrbit orbit;
         astro::JulianDate tt(orbit.dateFromSeconds(current_time));

         // Convert times to tdb since pulsar times are usually given in tdb time system
         double tdb = tt + orbit.tdb_minus_tt(tt)/86400.;
         double dt = (tdb - m_t0[i])*86400;

         // Calculate a dimensionless target number based on the Poisson distribution
         double target = -std::log(1.-RandFlat::shoot(1.));

         double area = 1.0e4 * EventSource::totalArea();
         double current_period = period(current_time,i);

         double num_per_cycle = m_flux[i]*area*current_period;
         double num_cycles = floor(target/num_per_cycle);

         m_interval[i] = current_period * num_cycles;

         // Determine what fraction of a cycle the pulsar is at for the current time
         double cycle_fraction = fmod(m_freq[i]*dt + 0.5 * m_freq_dot[i]*dt*dt,1);
             
         int phase_index = (int) floor(cycle_fraction * 21.);
         if(phase_index == 21) 
            phase_index = 0;

         double current = num_per_cycle * num_cycles;

         // Start sum by subtracting off part of bin that isn't used
         double phase_sum = - m_lightCurve[i][phase_index] * (cycle_fraction * 21. - 1.0 * phase_index) 
                            * area * (current_period / 21.);
         
         // Find out which bin corresponds to the target
         for(int j = 0; j < 43; j++)
         {
            // For each bin add rate * bin-time 
            phase_sum += m_lightCurve[i][phase_index] * area * (current_period / 21.);

            if(current + phase_sum > target)
            {
               m_interval[i] += (j+1.)*current_period / 21.;
               m_interval[i] -= current_period / 21. * (current + phase_sum - target) / (m_lightCurve[i][phase_index] * area * (current_period / 21.));
               break;
            }
            else
            {
               phase_index++; 
               if(phase_index == 21) phase_index = 0;
            }
         }

      }
      else
         m_interval[i] -= time_decrement;
   }
} 


double GalPulsars::period(double time, int pulsarIndex) const {
   return 1./m_freq[pulsarIndex] + (-m_freq_dot[pulsarIndex]/(m_freq[pulsarIndex]*m_freq[pulsarIndex]))*(time - m_t0[pulsarIndex]);
}

std::pair<double,double> GalPulsars::dir(double energy)
{
   if(m_changed == -1)
   {
      ;// Throw an exception
   }

   // return direction in galactic coordinates
   return std::make_pair<double,double>(m_lon[m_changed],m_lat[m_changed]);
}

float GalPulsars::operator()(float xi) const {

   if(m_changed == -1)
   {
      ; // Throw an exception
   }

    // single power law, or first segment
    return static_cast<float>(power_law(xi, m_lowCutoff[m_changed], m_highCutoff[m_changed], 1. + m_spectralIndex[m_changed]));
}

double GalPulsars::energy(double time) {
   return (*this)(RandFlat::shoot());
}


// differential rate: return energy distrbuted as e**-gamma between e1 and e2, if r is uniform from 0 to1
double GalPulsars::power_law( double r, double e1, double e2, double gamma) const
{
   return gamma==1
           ?  e1*exp(r*log(e2/e1))
           :  e1*exp(log(1.0 - r*(1.-pow(e2/e1,1-gamma)))/(1-gamma));
}

// verify that light curve agrees with flux for each source
void GalPulsars::initLightCurve(void)
{
   for(int i = 0; i < m_lightCurve.size(); i++)
   {
      double rate = 0;
      for(int j = 0; j < 21; j++)
         rate += m_lightCurve[i][j] / 21.;

      double scaling_factor = m_flux[i] / rate;

      for(int k = 0; k < 21; k++)
         rate *= scaling_factor;
   }
}