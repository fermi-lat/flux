// $Header: /nfs/slac/g/glast/ground/cvs/flux/src/test/testMgr.cxx,v 1.1.1.1 2003/07/29 18:22:19 burnett Exp $

//#include "FluxSvc/ISpectrumFactory.h"

#include "flux/EventSource.h"
#include "flux/SpectrumFactoryTable.h"
#include "flux/FluxMgr.h"

#include <iostream>
#include <fstream>
#include <algorithm>


static int default_count = 10 ;
//Testing
static const char * default_source="default";
//static const char * default_source="galdiffusemap";
//Default
//static const char * default_source="CrElectron";

void help() {
    std::cout << 
        "   Simple test program to create a particle source, then run it.\n"
        "   Command line args are \n"
        "      <source name> <count>\n"
        "   with defaults \"" 
        <<  default_source << "\"," << default_count
        << "\n  Also, 'help' for this help, 'list' for a list of sources and Spectrum objects "
        << std::endl;
}
void listSources(const std::list<std::string>& source_list ) {
    std::cout << "List of available sources:" << std::endl;
    for( std::list<std::string>::const_iterator it = source_list.begin(); 
        it != source_list.end(); ++it) { 
            std::cout << '\t'<< *it << std::endl;
        }

}
void listSpectra() {
    std::cout << "List of loaded Spectrum objects: " << std::endl;
    std::list<std::string> spectra(SpectrumFactoryTable::instance()->spectrumList());
    for( std::list<std::string>::const_iterator it = spectra.begin(); 
        it != spectra.end(); ++it) { 
            std::cout << '\t'<< *it << std::endl;
        }
}


#define DECLARE_SPECTRUM(x)   extern const ISpectrumFactory& x##Factory; x##Factory.addRef();

void flux_load() {

    // these are the spectra that we want to make available
    DECLARE_SPECTRUM( CHIMESpectrum);
    DECLARE_SPECTRUM( AlbedoPSpectrum);
    DECLARE_SPECTRUM( FILESpectrum);
    DECLARE_SPECTRUM( GalElSpectrum);
    DECLARE_SPECTRUM( SurfaceMuons);
    DECLARE_SPECTRUM( MapSpectrum);
 }

int main(int argn, char * argc[]) {
    using std::cout;
    using std::endl;
    flux_load();

    int count = default_count;
    std::string source_name(default_source);

    //TESTING MULTIPLE XML INPUT
    std::vector<std::string> fileList;
    fileList.push_back("$(FLUXROOT)/xml/source_library.xml");
    FluxMgr fm(fileList);

    //FluxMgr fm;

    //Testing the addfactory function
    //    static PencilBeam* sean=TestSpec::instance();
    //    fm.addFactory("seantest", sean );
    //End Test

    if ( argn >1 ) source_name = argc[1];
    if( source_name =="help") { help(); return 0; }
    if( source_name =="list") { 
        listSources(fm.sourceList());
        listSpectra(); return 0; }
    if ( argn >2 ) count = ::atoi(argc[2]);

    cout << "------------------------------------------------------" << endl;
    cout << " Flux test program: type 'help' for help" << endl;
    cout << ( ( argn ==1)?  " No command line args, using default fluxes \""
        :  " Selected source name \"");
    cout  << source_name <<"\"" << endl;


    std::list<std::string> source_list(fm.sourceList());

    if(( argn !=1) && std::find(source_list.begin(), source_list.end(), source_name)==source_list.end() ) {
        std::list<std::string> spectra(SpectrumFactoryTable::instance()->spectrumList());

        if( std::find(spectra.begin(), spectra.end(), source_name)==spectra.end() ) {
            std::cout << "Source \"" << source_name << "\" not found in the list or sources!" << std::endl;
            listSources(source_list);
            std::cout << "or in spectra list, which is:\n";
            listSpectra();

            return -1;
        }
    }
    std::list<std::string> allTheSources = fm.sourceList();
    std::list<std::string>::iterator abc;
    if(argn != 1){
        fm.test(std::cout, source_name, count);
        return 0;
    }
    std::string testfilename("testMgrOutput.out");
    std::ostream* m_out = new std::ofstream(testfilename.c_str());
    std::ostream& out = *m_out;
    std::cout << "Writing test results to the file " << testfilename << std::endl;
    
    for(abc= allTheSources.begin() ; abc != allTheSources.end() ; abc++){

        // now have FluxMgr create and run it.
        std::cout << "Source:  " << *abc << std::endl;
        out << "Source:  " << *abc <<std::endl;
        fm.test(out, (*abc), count);
        out << std::endl << std::endl << std::endl;
    }

    return 0;    
}

