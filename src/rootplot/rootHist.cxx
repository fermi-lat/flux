/** @file rootHist.cxx

$Header$
*/
//  class rootHist
//  Author:  Theodore Hierath
//  This class contains a basic histogram function
//  The bins are numbered from 0 to number of bins - 1;

#include <iostream>
#include "rootHist.h"
#include "flux/FluxException.h"

// Specify the number of bins for the histogram
rootHist::rootHist(int bins) : num_bins(bins)
{
    if(bins <= 0)
    {
        FATAL_MACRO( "ERROR in constructor for roothist:\n" 
            "   The total number of bins must be greater than zero.");
    }
    
    hist = new double[num_bins];
    for (int i = 0; i < num_bins; i++) hist[i] = 0;
}

rootHist::~rootHist(void)
{
    delete [] hist;
    hist = NULL;
}

rootHist::rootHist(const rootHist &oldHist) : num_bins(oldHist.num_bins)
{
    hist = new double[num_bins];
    for (int i = 0; i < num_bins; i++) hist[i] = oldHist.hist[i];
}

// Update a specific bin by replacing its contents
void rootHist::updateBin(int binNumber, double data)
{
    if(binNumber < 0) binNumber=0;
    else if(binNumber>=num_bins) binNumber= num_bins-1;
    hist[binNumber] = data;
}

// Retrieve the contents of a bin
double rootHist::retrieveBin(int binNumber)
{
    if(binNumber < 0) binNumber=0;
    else if(binNumber>=num_bins) binNumber= num_bins-1;
    return hist[binNumber];
}

// Increment the contents of a bin
void rootHist::incrementBin(int binNumber)
{
   if(binNumber < 0) binNumber=0;
   else if(binNumber >= num_bins) binNumber = num_bins-1;
   hist[binNumber]++;
}