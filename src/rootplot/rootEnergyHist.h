//  class rootEnergyHist
//  Author:  Theodore Hierath
//  This class contains a histogram and graphing class for
//  energies.
//  The bins are numbered from 0 to number of bins - 1;

#ifndef ROOT_ENERGY_HIST
#define ROOT_ENERGY_HIST

#include "rootHist.h"
#include <cmath>
#include <string>
#include <iostream>

/** This class stores and graphs the flux vs.\ energy.
    Normal usage typically entails
    1. calling the constructor
    2. calling setTitle, setXLabel, setYLabel, etc.
    3. storing energies in the object using store
    4. calling apply to normalize the histogram
    5. calling draw or writeFile to output the results
    @author Theodore Hierath
  */

class rootEnergyHist 
{
public:
   /** Constructor creates an empty energy histogram
      @params bins The number of bins to use in the graph
      @params min_energy The minimum energy in GeV for the binning range
      @params max_energy The maximum energy in GeV for the binning range   */
   rootEnergyHist(int bins, double min_energy, double max_energy);
   ~rootEnergyHist(void);
   rootEnergyHist(const rootEnergyHist& oldHist);

   /// Sets the title for the graph
   void setTitle(std::string title);

   /// Sets the x-axis label on the flux vs. energy graph
   void setXLabel(std::string label);

   /// Sets the y-axis label on the flux vs. energy graph
   void setYLabel(std::string label);

   /** Sets the type of scaling used on the graph
      Valid values are "linear", "semilogx", "semilogy", and "log"    */
   void setGraphType(const char *graph_type);

   /// Increments the temporary bin (by +1) corresponding to the given energy
   void store(double energy);

   /// Retrieves the raw count associated with the given bin number
   double retrieveCount(int binNumber);

   /// Retrieves the flux associated with the given bin number
   double retrieveFlux(int binNumber);

   /// Retrieves the flux uncertainty associated with the given bin number
   double retrieveFluxUncertainty(int binNumber);

   /// Retrieves the energy associated with the given bin number
   double retrieveEnergy(int binNumber);

   /// Retrieves the logorithmic (base 10) energy range for the histogram
   double retrieveRange(void);

   /// This function changes the graph type from E*flux vs. E to flux vs. E
   void setFluxMode(void);

   /// This function manually sets the minimum flux to be displayed on the graph
   void setFluxMin(double f);

   /// This function manually sets the maximum flux to be displayed on the graph
   void setFluxMax(double f);

   /** Applies the scaling factor to the raw counts to obtain the fluxes.
       This function clears the raw count histogram, stores the flux values
       in the flux histogram, and calculates the uncertainties.  If there was
       a previous value in the flux and flux uncertainty histograms, this function
       does not overwrite the values, but adds them together.  This function must
       be used before the draw member function.   */
   void apply(double scale_factor);

   /// This function clears everything in the object except for the number of bins.
   void reset(void);

   /** Writes the contents of the flux and flux uncertainty histograms to a script file.
       This file is then interpreted by root.
       @param scale_factor This is a final scale factor applied to the graphs.
       @param mode This can be "normal", "begin", or "end".  Normal means that the
          draw function for this object will be used alone.  Begin means that this
          draw function is first called before the corresponding one for rootEnergyHist.
          End means that this draw function is first called after the draw function
          for rootEnergyHist.
       @param current_plot This is the number of times this function has been called.
          If it is the first time it is called, this number should be '0'.  This is 
          used by the function so that it knows if it is at the beginning, middle, or
          end of the script file.
       @param total_plots  This is the total number of plots that will be graphed.  This
          is used in conjunction with current_plot so that the function knows when the 
          script file should end.   */
   void draw(double scale_factor, std::string mode, int current_plot, int total_plots);

   /** Writes the contents of the flux histogram to an output file with the following
       format.  The first line is a comment describing the file.  The following lines
       contain two pairs of numbers, the energy (MeV), and the flux 
       (particles/(m^2*s*sr*MeV)).
       @param scale_factor The final scale factor to apply
       @param out_file The output file stream to write to   */
   void writeFile(double scale_factor, std::ostream& out_file);
private:
   rootHist *flux_hist; ///< histogram to store fluxes
   rootHist *flux_sigma_hist; ///< histogram to store uncertainties
   rootHist *raw_hist; ///< histogram to store raw counts
   const int num_bins;  ///< number of bins for the histogram
   const double emin;   ///< minimum energy
   const double emax;   ///< maximum energy
   double range;        ///< logarithmic range of energies
   std::string graphTitle;    ///< title of graph
   std::string xlabel;  ///< x axis title
   std::string ylabel;  ///< y axis title 

   /// Enumeration that determines the graph type
   enum graphType
   {
      linear,     ///< Linear graph
      semilogx,   ///< Linear y-axis, logarithmic x-axis
      semilogy,   ///< Linear x-axis, logarithmic y-axis
      loglog      ///< Logarithmic x and y axes
   } currentType;

   bool fluxMode;  ///< If true flux is graphed instead of E*flux
   bool use_flux_min; ///< If true flux_min is used as the minimum flux on the graph
   bool use_flux_max; ///< If true flux_max is used as the maximum flux on the graph
   double flux_min; ///< Manually specified minimum flux
   double flux_max; ///< Manually specified maximum flux
};

#endif // ROOT_ENERGY_HIST
