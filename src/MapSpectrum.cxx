
#include "flux/MapSpectrum.h"

#include <cmath>
#include <algorithm>
#include <functional>
#include <fstream>
#include <iostream>

#include "fitsio.h"  
#include "CLHEP/Random/RandFlat.h"

// this is needed to include in the executable or dll
#include "flux/SpectrumFactory.h"

static SpectrumFactory<MapSpectrum> factory;
const ISpectrumFactory& MapSpectrumFactory = factory;

namespace {

// This routine breaks down a string into its components based on the
// characters appearing in the delimiters string.  A vector
// of strings is returned.
   std::vector<std::string> string_split(std::string input, 
                                         const std::string &delimiters) {
      std::vector<std::string> components;
      
      std::string::size_type j;
      while ( (j = input.find_first_of(delimiters)) != std::string::npos ) { 
         if (j != 0) components.push_back(input.substr(0, j));
         input = input.substr(j+1);
      }
      components.push_back(input);
      return components;
   }
} // unnamed namespace

MapSpectrum::MapSpectrum()
:m_mapE0(10.),m_extrapE0(10.),m_EMax(1000.)
{}

MapSpectrum::MapSpectrum(const std::string& params)
:m_mapE0(parseParamList(params,0)),m_extrapE0(parseParamList(params,1)),
	m_EMax(parseParamList(params,2)),m_defaultIndex(parseParamList(params,3))
{
    m_binSize = 0.5; //this should come from FITS file instead!
    //let's set the map to zero right off the bat:
    for(double l=-180 ; l<=180 ; l+=m_binSize){
        for(double b=-90 ; b<=90 ; b+=m_binSize){
			cellinfo temp;
			temp.intensity = 0;
		    m_testCatalog.push_back(temp);
       }
    }

    //also initialize the row-counting map:
    for(double b=-90 ; b<=90 ; b+=m_binSize){
            m_rowCatalog[RowKey(b,m_binSize)] = 0;
        }

    m_particle_name = "gamma";
    
// Parse the parameter list given in the params string.
//
// Put each element of params into a vector using "," as the
// delimiter.
   std::vector<std::string> paramVector = ::string_split(params, ",");

// Assume the fifth element, if it exists, is the map filename,
// otherwise use the default.
   if (paramVector.size() == 5) {
      initialization_document = paramVector[4];
   } else {
      initialization_document = "src/sources/gas.gal.gz";
   }
    
    int status=0;
    int fitsflag = 1; //1 if we have a FITS file, 0 if not
    const char* flux_root = ::getenv("FLUXROOT");
    //define the file
    std::string doc_path2= (flux_root? std::string(flux_root)+"/" : "");
    std::string fileName2 = doc_path2+initialization_document;
    //----------------------fitsio----------------------------------------------------
    //open the file
    fitsfile *fptr;
    ffopen(&fptr, fileName2.c_str(), 0, &status);
    //printf("\nOpening file: %s status = %d\n", fileName2.c_str(),status);
    if(!status){
        //   fitsflag = 0; //this is not a proper FITS file
        
        //ok, the FITS file has been opened.  Now, get the data:
        int nfound, anynull;
        long naxes[2], fpixel, nbuffer, npixels, ii;
        long jj=0;
        
#define buffsize 720
        float datamin, datamax, nullval;
        float buffer[buffsize];

        int buffersize = 360*m_binSize;

        status = 0;
        
        /* read the NAXIS1 and NAXIS2 keyword to get image size */
        if ( fits_read_keys_lng(fptr, "NAXIS", 1, 2, naxes, &nfound, &status) )
            std::cout << "error: " << status << std::endl;
        //std::cout << "status: " << status << std::endl;

        m_binSize = 360./naxes[0];

        npixels  = naxes[0] * naxes[1];         /* number of pixels in the image */
        fpixel   = 1;
        nullval  = 0;                /* don't check for null values in the image */
        datamin  = (float)1.0E30;
        datamax  = (float)-1.0E30;
        /*double */float curl,curb,curint,curind;
        curl = -180.-m_binSize;
        curb = -90.;//-m_binSize;
        while (npixels > 0)
        {

            nbuffer = npixels;
            if (npixels > buffsize)
                nbuffer = buffsize;     /* read as many pixels as will fit in buffer */
            
            /* Note that even though the FITS images contains unsigned integer */
            /* pixel values (or more accurately, signed integer pixels with    */
            /* a bias of 32768),  this routine is reading the values into a    */
            /* float array.   Cfitsio automatically performs the datatype      */
            /* conversion in cases like this.                                  */
            
            if ( fits_read_img(fptr, TFLOAT, fpixel, nbuffer, &nullval,
                buffer, &anynull, &status) )
                std::cout << "error in reading image: " << status << std::endl;
            
			//initialize the key structure:
			Key abc(0,0,0.5);
            for (ii = 0; ii < nbuffer; ii++)  {
                //each execution of this loop should move one bin "right"
                curl += m_binSize;
                if(curl == 180){
                    curl = -180;
                    curb += m_binSize;
                }
                curint=buffer[ii];
                curind = m_defaultIndex;

				abc.init(curl,curb,m_binSize);
			    int trig=abc.operator int();
				//the 1E4 is due to the map being held in 
				//particles/s*cm^2*sr in stead of p/s*m^2*sr
				//do the extrapolation down to lower energy here
			    m_testCatalog[trig].intensity = curint*(pow(double(m_extrapE0/m_mapE0),(1.0-curind)))*1E4;
				m_testCatalog[trig].index = curind;
          
                //          std::cout << (int)curl << ", " << (int)curb << ", " << m_catalog[std::make_pair<int,int>(curl,curb)].intensity*100000. << std::endl;
                
            }
            jj++;
            npixels -= nbuffer;    /* increment remaining number of pixels */
            fpixel  += nbuffer;    /* next pixel to be read in image */
        }
        //------------------------end fitsio input---------------------------------
    }else{
        //-----------------------begin ascii input----------------------------------
        //ok, this is not a proper FITS file.  Maybe it's in the Ascii format?
        std::ifstream input_file;
        input_file.open(fileName2.c_str());
        
        if(false == input_file.is_open())
        {
            std::cerr << "ERROR:  Unable to open:  " << fileName2.c_str() << std::endl;
            exit(0);
        }
        else
        {
            double curl,curb,curint,curind;
            //initialize the key structure:
			Key abc(0,0,0.5);
            while (!input_file.eof()){
                input_file >> curl;
                input_file >> curb;
                input_file >> curint;
                input_file >> curind;

				abc.init(curl,curb,m_binSize);
			    int trig=abc.operator int();

				//the 1E4 is due to the map being held in 
				//particles/s*cm^2*sr in stead of p/s*m^2*sr
			    m_testCatalog[trig].intensity = curint*(pow(double(m_extrapE0/m_mapE0),(1.0-curind)))*1E4;
				m_testCatalog[trig].index = curind;
            }
        }
    }
    //-------------------------end ascii input-----------------------------------
    
    //either way, now figure out the total flux.
    setNetFlux();

// Assume the fourth element of paramVector, if it exists, is the
// net flux to be used.

   if (paramVector.size() == 4) {
      float new_netFlux = ::atof(paramVector[3].c_str());

// Since m_rowCatalog is used as the cumulative distribution for
// drawing the direction of the particle, we must rescale to the new
// value of the net flux.  This must be done for m_testCatalog as
// well.
      
      for (std::map<RowKey, double>::iterator it = m_rowCatalog.begin();
           it != m_rowCatalog.end(); it++) {
         it->second *= new_netFlux/m_netFlux;
      }
      for (std::vector<cellinfo>::iterator it = m_testCatalog.begin();
           it != m_testCatalog.end(); it++) {
         it->intensity *= new_netFlux/m_netFlux;
      }
      m_netFlux = new_netFlux;
   }

    //std::cout << "net flux is " << m_netFlux << std::endl;
    if ( fits_close_file(fptr, &status) )
        std::cout << "error in closing file: " << status << std::endl;
    
    }
    
    /// calculate flux for the current position
    double MapSpectrum::flux(double time) const
    {
        return m_netFlux;
    }
    
    std::pair<float,float> MapSpectrum::dir(float energy)const{
        return std::make_pair<float,float>(0.0,0.0);
    }
    
    double MapSpectrum::energy( double time){
		findDir(0.); //energy here is just a dummy - now we have direction, index.
        return operator()(RandFlat::shoot());
    }
    
    
    
    ///this returns a galactic direction, in the form (l,b)
    std::pair<double,double> MapSpectrum::findDir(double energy){
        //here is where we decide where the next particle will come from, and set the associated flux and energy of the particle.
        double remainingFlux=RandFlat::shoot(m_netFlux);

        double l=-180;
        double b=-90;
        
		//ok, now count over the rows until there isn't enough remaining flux
		//to account for a whole row.  This is then the row to use:
        while(remainingFlux >= m_rowCatalog[RowKey(b,m_binSize)]){
            remainingFlux -= m_rowCatalog[RowKey(b,m_binSize)];
            b+=m_binSize;
        }

		//initialize the key structure:
		Key abc(0,0,0.5);
		Key abc2(0,0,0.5);

        //get the l,b coordinates:  count over bins until no more flux remains
        while(remainingFlux /*>=*/ > 0){

			abc.init(l,b,m_binSize);
			int trig=abc.operator int();
			remainingFlux -= m_testCatalog[trig].intensity*sizeOfBin(b);

            l+=m_binSize;

        }
        //we have gone one bin too far, in the counting - need to back up one:
		l -= m_binSize;

        //now set the flux:
        m_flux = m_netFlux;
        
        //and the energy index of this particle:
		abc2.init(l,b,m_binSize);
		int trig2=abc2.operator int();
		m_index = m_testCatalog[trig2].index;

		//now that we've gotten the bin from the map, 
		//make sure the photon comes from somewhere in that bin:
		double spreadl = (RandFlat::shoot()-0.5)*m_binSize;
		double spreadb = (RandFlat::shoot()-0.5)*m_binSize;

		//now set the current directions:
		m_currentl=l+spreadl;
		m_currentb=b+spreadb;

        return std::make_pair<double,double>(l,b);    
    }
    
	///this returns a galactic direction, in the form (l,b)
    std::pair<double,double> MapSpectrum::dir(double energy){
		//we should have already made these, just call them:
		return std::make_pair<double,double>(m_currentl,m_currentb);
	}
    
    /// sample a single particle energy from the spectrum
    float MapSpectrum::operator() (float f)const
    {
		double index = m_index;
        //double m_emax = 100;
        if( index == 0.0 )     return m_extrapE0;
        
        if( index == 1.0 ) return m_extrapE0*exp(f*log(m_EMax/m_extrapE0));
        
        float x = 1 - exp((1-index)*log(m_EMax/m_extrapE0));
        return m_extrapE0*exp(log(1-x*f)/(1-index)); //the 1000 is for Mev/Gev conversion
    }
    
    
    std::string MapSpectrum::title() const
    {
        return "MapSpectrum";
    }
    
    const char * MapSpectrum::particleName() const
    {
        return m_particle_name.c_str();
    }
    
	/// This function sets the cumulative flux over the whole sky, and 
	/// also in each row, for later integration.
    void MapSpectrum::setNetFlux(){
        m_netFlux=0.;
        double cumFlux = 0.;//the cumulative flux, to be used as an index.
		//initialize the key structure:
		Key abc(0,0,0.5);
        for(double l=-180 ; l<180 ; l+=m_binSize){        
            //double rowFlux = 0.;//the flux summed over the current row (constant l)
            for(double b=-90 ; b<90 ; b+=m_binSize){
	
				abc.init(l,b,m_binSize);
		        int trig=abc.operator int();
		        m_index = m_testCatalog[trig].index;

				cumFlux+=m_testCatalog[trig].intensity*sizeOfBin(b);
                m_netFlux+=m_testCatalog[trig].intensity*sizeOfBin(b);
                m_rowCatalog[RowKey(b,m_binSize)] += m_testCatalog[trig].intensity*sizeOfBin(b);
                            
			}
        }
        return;
    }
    
    float MapSpectrum::parseParamList(std::string input, int index)
    {
        std::vector<float> output;
        int i=0;
        for(;!input.empty() && i!=std::string::npos;){
            float f = ::atof( input.c_str() );
            output.push_back(f);
            i=input.find_first_of(",");
            input= input.substr(i+1);
        } 
        return output[index];
    }
    
    double MapSpectrum::sizeOfBin(double b){
        double binTop = (b+(m_binSize/2.))*M_PI/180.;
        double binBot = (b-(m_binSize/2.))*M_PI/180.;

        // integrated solid angle of bin(in steradians):
        return m_binSize*(sin(binTop)-sin(binBot))*(M_PI/180.);

    }
    
    double MapSpectrum::interval(double time){return -1;}
