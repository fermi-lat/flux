//  class rootHist
//  Author:  Theodore Hierath
//  This class contains a basic histogram function
//  The bins are numbered from 0 to number of bins - 1;

#ifndef ROOT_HIST
#define ROOT_HIST

/** This class serves as an array that checks for invalid indexing.
    @author Theodore Hierath
*/

class rootHist
{
public:
   /// Constructor takes number of bins for histogram.
   rootHist(int bins);
   rootHist(const rootHist& oldHist);
   ~rootHist(void);

   /**  Stores a number in the specified bin
        @params binNumber The bin to store the data in.  The bin numbering
                goes from 0 to (total number of bins - 1).
        @params data The number to store   */
   void updateBin(int binNumber, double data);

   /** Returns the contents of the specified bin.  The bin numbering
       goes from 0 to (total number of bins - 1).  */
   double retrieveBin(int binNumber);

   /** Increments the value stored in the bin by 1.
       @params binNumber The bin to increment.  The bin numbering goes
               from 0 to (total number of bins - 1). */
   void incrementBin(int binNumber);

private:
   double *hist;  ///< The histogram
   const int num_bins;  ///< The number of bins in the histogram
};

#endif // ROOT_HIST