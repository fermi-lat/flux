//  class rootAngleHist
//  Author:  Theodore Hierath
//  This class contains a histogram and graphing class for
//  angles.
//  The bins are numbered from 0 to number of bins - 1;

#ifndef ROOT_ANGLE_HIST
#define ROOT_ANGLE_HIST

#include "rootHist.h"
#include <string>
#include <iostream>

/**
   This class stores the the flux vs.\ incident angle information
   for creating a graph.
   Normal usage typically entails
    1. calling the constructor
    2. calling setTitle, setXLabel, setYLabel, etc.
    3. storing energies in the object using store
    4. calling apply to normalize the histogram
    5. calling draw or writeFile to output the results
   @author Theodore Hierath
*/   
class rootAngleHist 
{
public:
   /// Constructor takes number of bins for histogram.
   rootAngleHist(int bins);
   ~rootAngleHist(void);
   rootAngleHist(const rootAngleHist& oldHist);

   /// Sets the title for the graphs
   void setTitle(std::string title);

   /// Sets the x-axis label on the flux vs. azimuth angle graph
   void setPhiXLabel(std::string label);

   /// Sets the y-axis label on the flux vs. azimuth angle graph
   void setPhiYLabel(std::string label);

   /// Sets the x-axis label on the flux vs. zenith angle graph
   void setThetaXLabel(std::string label);

   /// Sets the y-axis label on the flux vs. zenith angle graph
   void setThetaYLabel(std::string label);

   /** Sets the y-axis scaling for the graphs
       Valid values are "linear" and "log"   */
   void setGraphType(const char *graph_type);
   
   /// Increments the temporary bin corresponding to parameter theta
   void storeTheta(double theta_value);

   /// Increments the temporary bin corresponding to parameter phi
   void storePhi(double phi_value);

   /// Retrieves the phi angle raw count from the specified bin
   double retrievePhiCount(int binNumber);

   /// Retrieves the phi angle flux from the specified bin
   double retrievePhiFlux(int binNumber);

   /// Retrieves the phi angle flux uncertainty from the specified bin
   double retrievePhiFluxUncertainty(int binNumber);

   /// Retrieves the theta angle raw count from the specified bin
   double retrieveThetaCount(int binNumber);

   /// Retrieves the theta angle flux from the specified bin
   double retrieveThetaFlux(int binNumber);

   /// Retrieves the theta angle flux uncertainty from the specified bin
   double retrieveThetaFluxUncertainty(int binNumber);

   /// Retrieves the angle in degrees associated with the bin number
   double retrievePhiAngle(int binNumber);

   /// Retrives the value of cos(theta) associated with the bin number
   double retrieveThetaAngle(int binNumber);

   /** Applies the scaling factor to the raw counts to obtain the fluxes.
       This function clears the raw count histograms, stores the flux values
       in the flux histograms, and calculates the uncertainties.  If there was
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

   /// Not yet implemented
   void writeFile(double scale_factor, std::ostream& out_file);
private:
   rootHist *theta_hist; ///< histogram for zenith angle fluxes
   rootHist *theta_sigma_hist;  ///< histogram to store zenith angle uncertainties
   rootHist *theta_raw_hist;  ///< histogram to store zenith agnle raw counts
   rootHist *phi_hist;   ///< histogram for azimuthal angle fluxes
   rootHist *phi_sigma_hist;  ///< histogram to store azimuthal angle uncertainties
   rootHist *phi_raw_hist;  ///<  histogram to store azimuthal angle raw counts
   const int num_bins;  ///< number of bins for the histogram
   std::string graphTitle;    ///< title of graph
   std::string PhiXLabel;  ///< phi graph x axis title
   std::string PhiYLabel;  ///< phi graph y axis title 
   std::string ThetaXLabel; ///< theta graph x axis title
   std::string ThetaYLabel; ///< theta graph y axis title
   enum graphType{linear,log} currentType; ///< Determines whether the y-axis will be linear or logarithmic
};

#endif // ROOT_ANGLE_HIST
