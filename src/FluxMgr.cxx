/** @file FluxMgr.cxx
    @brief Implementation of FluxMgr

  $Header: /nfs/slac/g/glast/ground/cvs/flux/src/FluxMgr.cxx,v 1.1.1.1 2003/07/29 18:22:19 burnett Exp $
*/

#include "flux/FluxMgr.h"
#include "flux/EventSource.h"
#include "flux/SpectrumFactoryTable.h"
#include "flux/GPS.h"
#include "flux/FluxException.h" // defines FATAL_MACRO
#include "CompositeSource.h"

#include <xercesc/dom/DOM_Document.hpp>
#include <xercesc/dom/DOM_Element.hpp>
#include "xml/Dom.h"
#include "xml/IFile.h"

#include <sstream>

#define DLL_DECL_SPECTRUM(x)   extern const ISpectrumFactory& x##Factory; x##Factory.addRef();
#define DECLARE_SPECTRUM(x)   extern const ISpectrumFactory& x##Factory; x##Factory.addRef();



FluxMgr::FluxMgr(const std::vector<std::string>& fileList, std::string dtdname)
: m_dtd(dtdname.empty()? "$(FLUXROOT)/xml/source.dtd" : dtdname )
{
    if( fileList.empty() ){
        defaultFile();
    }else{		
        init(fileList);		
    }	
}

void FluxMgr::defaultFile(){
    //Purpose: to set the default xml file and initialize the package to use it.
    std::vector<std::string> input;

    // must find the source_library.xml file.
    // set up the xml document to use for initialization parameters
    const char* flux_root = ::getenv("FLUXROOT");
    std::string doc_path= (flux_root? std::string(flux_root)+"/xml/" : "");
    input.push_back(/*doc_path+initialization_document*/"$(FLUXROOT)/xml/source_library.xml");	
    init(input);
}

void FluxMgr::init(const std::vector<std::string>& fileList){	
    
    //Purpose:  to initialize the infrastructure of the package with
    //the named xml files.
    std::string fileName;
    
    xml::XmlParser parser;
    
    std::string xmlFileIn = writeXmlFile(fileList);
    
    // a quick way of displaying what goes to the parser
    //std::cout << xmlFileIn <<std::endl;
    
    m_library_doc = parser.parse(xmlFileIn);
    
    if (m_library_doc == DOM_Document()) {
        FATAL_MACRO("Parse error: processing the document" << std::endl
            << xmlFileIn << std::endl);
        return;
    }
    
    // Root element is of type source_library.  Content is
    // one or more source elements.
    
    s_library = m_library_doc.getDocumentElement();
    
    // loop through the source elements to create a map of names, DOM_Elements
    if (s_library != DOM_Element()) {
        
        DOM_Element child = xml::Dom::getFirstChildElement(s_library);
        DOM_Element toplevel = xml::Dom::getFirstChildElement(s_library);
        
        while (child != DOM_Element()) {
            while (child.getAttribute("name") == DOMString()) {
                s_library = child;
                child = xml::Dom::getFirstChildElement(s_library);
            }
            
            while (child != DOM_Element()) {
                std::string name = xml::Dom::transToChar(child.getAttribute("name"));
                //std::cout << name << std::endl;
                std::string parentfilename = xml::Dom::transToChar(toplevel.getAttribute("title"));
                m_sources[name]=std::make_pair<DOM_Element,std::string>(child,parentfilename);
                child = xml::Dom::getSiblingElement(child);
            }
            
            child = xml::Dom::getSiblingElement(toplevel);
            toplevel=child;
        }
        
    }
    // these are the locally defined spectra that we want to make available
    DLL_DECL_SPECTRUM( FILESpectrum);
    DLL_DECL_SPECTRUM( TimeCandle);

    DLL_DECL_SPECTRUM( SurfaceMuons);

    // these are deprecated, will be replaced by Hiroshima group
    DLL_DECL_SPECTRUM( AlbedoPSpectrum);
    DLL_DECL_SPECTRUM( CHIMESpectrum );
    DLL_DECL_SPECTRUM( GalElSpectrum);
    
}


FluxMgr::~FluxMgr(){
}

EventSource* FluxMgr::source(std::string name)
{
    //Purpose: to return a pointer to a source, referenced by name.
    //Input: the name of the desired source.
    // first check that it is in the library
    if( m_sources.find(name)==m_sources.end() ) {
        return 0;
    }
    return getSourceFromXML(m_sources[name].first);
}


EventSource*  FluxMgr::getSourceFromXML(const DOM_Element& src)
{
    //Purpose: sourceFromXML - create a new EventSource from a DOM element
    //instantiated, e.g., from a description in source_library.xml
    //Input:  the element holding particle information.
    DOM_Node    childNode = src.getFirstChild();
    if (childNode == DOM_Node()) {
    /*
    FATAL_MACRO("Improperly formed XML event source");
    return 0;
        */
        // no child node: expect to find the name defined.
        return  new FluxSource(src);
    }
    
    DOM_Element sname = xml::Dom::getFirstChildElement(src);
    if (sname == DOM_Element() ) {
        FATAL_MACRO("Improperly formed XML event source");
        return 0;
    }
    // If we got here, should have legit child element
    if ((sname.getTagName()).equals("spectrum")) {
        return  new FluxSource(src);
    }
    else if ((sname.getTagName()).equals("nestedSource")) {
        
        // Search for and process immediate child elements.  All must
        // be of type "nestedSource".  There may be more than one.
        // Content model for nestedSource is EMPTY, so can omit check
        // for that in the code
        
        CompositeSource* cs;
            cs = new CompositeSource();
        do { 
            DOM_Element selem = 
                getLibrarySource(sname.getAttribute("sourceRef"));
            if (selem == DOM_Element()) {
                FATAL_MACRO("source name" << 
                    xml::Dom::transToChar(sname.getAttribute("sourceRef")) << 
                    "' not in source library");
            }
            cs->addSource(getSourceFromXML(selem)); 
            sname = xml::Dom::getSiblingElement(sname);
        } 
        while (sname != DOM_Element() );
        return cs;
    }
    else {
        FATAL_MACRO("Unexpected element: "<< 
            xml::Dom::transToChar(sname.getTagName()) );
    }
    return 0;
}




DOM_Element    FluxMgr::getLibrarySource(const DOMString& id)
{
    //Purpose: source library lookup.  Each source is uniquely identified
    // by its "name" attribute because "name" is of type ID
    
    // quit if the library was unitialized
    if (s_library == DOM_Element() ) return DOM_Element(); 
    
    return m_library_doc.getElementById(id);
}

std::list<std::string> FluxMgr::sourceList() const
{
    std::list<std::string> s;
    for( std::map<std::string, std::pair<DOM_Element,std::string> >::const_iterator it = m_sources.begin();
    it != m_sources.end();
    ++it){
        s.push_back((*it).first);
    }
    return s;
}

std::vector<std::pair< std::string ,std::list<std::string> > > FluxMgr::sourceOriginList() const
{
    std::vector<std::pair< std::string ,std::list<std::string> > > originList;
    for( std::map<std::string, std::pair<DOM_Element,std::string> >::const_iterator it = m_sources.begin();
    it != m_sources.end();
    ++it){
        //now see if we can find the filename already used in the vector:
        std::vector<std::pair< std::string ,std::list<std::string> > >::iterator topiter;
        for(topiter = originList.begin() ; topiter != originList.end() ; topiter++){
            if( (*topiter).first==(*it).second.second) break;
        }

        if( topiter != originList.end() ){ (*topiter).second.push_back((*it).first);
        }else{
            std::list<std::string> abc;
            abc.push_back((*it).first);
            originList.push_back(std::make_pair< std::string ,std::list<std::string> >((*it).second.second,abc) );
        }

    }

    return originList;
}

/// generate some test output
void FluxMgr::test(std::ostream& cout, std::string source_name, int count)
{   
    EventSource* e = source(source_name);
    setExpansion(1.);
    double time=0.;
    
    const int howMany = e->howManySources();
    //int counts[howMany] = 0;
    std::vector<int> counts;
    
    //std::vector<int>::const_iterator countIter = counts.begin();
    
    int q;
    for(q=0 ; q<=howMany+2 ; q++){
        counts.push_back(0);
        //  countIter++;
    }
    
    cout << "running source: " << e->fullTitle() << std::endl;
    cout << " Total rate is: " << e->rate(time) << " Hz into " << e->totalArea() << " m^2" << std::endl;
    //cout << "LaunchType" << f->retLaunch() << "Pointtype" << f->retPoint() <<std::endl;
    cout << " there are " << howMany << " Sources total..." << std::endl;
    cout << "    Generating " << count << " trials " << std::endl;
    cout << " --------------------------------" << std::endl;
    
    //testing rotateangles function
    GPS::instance()->rotateAngles(std::make_pair<double,double>(0.0,0.3));
    EventSource* f;
    double totalinterval=0;
    for( int i = 0; i< count; ++i) {
        
        f = e->event(time);
        //TESTING THE lat, lon FUNCTIONS
        //cout << std::endl << "lat=" << GPS::instance()->lat() << ' ' <<"lon=" << GPS::instance()->lon() << std::endl;
        //double curTime=GPS::instance()->time();
        //cout << std::endl << "testlat=" << GPS::instance()->orbit()->testLatitude(curTime) << ' ' << "testlon=" << GPS::instance()->orbit()->testLongitude(curTime) << std::endl;
        
        double interval=e->interval(time);
        
        //here we increment the "elapsed" time and the "orbital" time,
        //just as is done in flux.  NOTE: this is important for the operation 
        //of fluxsource, and is expected.
        time+=interval;
        pass(interval);
        int sourceNumber = e->numSource();
        if(sourceNumber==-1){counts[0]++;
        }else{counts[sourceNumber]++;}
        
        totalinterval+=interval;
        cout << f->particleName()
            << "(" << f->energy()<< " GeV)"
            << ", Launch: "  << f->launchPoint() 
            << ", Dir "      << f->launchDir() 
            << ", Flux="     << f->flux(time) 
            << ", Interval=" << interval ;
        if(sourceNumber!=-1) cout <<", SourceID: "<< sourceNumber ;
        cout << "\nElapsed time= " << totalinterval 
            << std::endl;
    }
    cout << "------------------------------------------------------" << std::endl;
    
    cout << std::endl << "Average Interval=" << totalinterval/count <<" , "
        << "Average rate = " << count/totalinterval <<std::endl;
    
    cout << "Source Statistics: " << std::endl;
    for(q=0 ; q<howMany ; q++){
        cout << "source #" << q+1 << ": " << counts[q] << " events counted." << std::endl;
    }
    
    
}


void FluxMgr::addFactory(std::string name, const ISpectrumFactory* factory ) {
    SpectrumFactoryTable::instance()->addFactory(name,factory);
}

void FluxMgr::setExplicitRockingAngles(std::pair<double,double> ang){
    GPS::instance()->rotateAngles(ang);
}

std::pair<double,double> FluxMgr::getExplicitRockingAngles(){
    return GPS::instance()->rotateAngles();
}

/// set the desired pointing history file to use:
void FluxMgr::setPointingHistoryFile(std::string fileName){
	GPS::instance()->setPointingHistoryFile(fileName);
}

void FluxMgr::setExpansion (double p){
    // set the expansion factor for the orbit (-1) = random
    GPS::instance()->expansion(p);
}

// pass a specific amount of time
void FluxMgr::pass(double t){
    GPS::instance()->pass(t);
    synch();
}

GPStime FluxMgr::time () const{
    return GPS::instance()->time();
}

void FluxMgr::synch(){
    GPS::instance()->synch();
}

void sampleintvl ( /*GPStime*/double t ){
    GPS::instance()->sampleintvl(t);
}

//get the current satellite location
std::pair<double,double> FluxMgr::location(){
    return std::make_pair<double,double>(GPS::instance()->lat(),GPS::instance()->lon());
}

//get the transformation matrix.
HepRotation FluxMgr::CELTransform(double time){
    return GPS::instance()->CELTransform(time);
}

//get the transformation matrix.
HepRotation FluxMgr::orientTransform(double time){
    return GPS::instance()->rockingAngleTransform(time);
}

///this transforms glast-local (cartesian) vectors into galactic (cartesian) vectors
HepRotation FluxMgr::transformGlastToGalactic(double time){
    return GPS::instance()->transformGlastToGalactic(time);
}

///this sets the rocking mode in GPS.
std::vector<double> FluxMgr::setRockType(GPS::RockType rockType, double rockAngle){
   int type=GPS::instance()->setRockType(rockType);
   double degrees = GPS::instance()->rockingDegrees(rockAngle);
   std::vector<double> ret;
   ret.push_back(type);
   ret.push_back(degrees);
   return ret;
}

///this sets the rocking mode in GPS.
std::vector<double> FluxMgr::setRockType(int rockType, double rockAngle){
   int type=GPS::instance()->setRockType(rockType);
   double degrees = GPS::instance()->rockingDegrees(rockAngle);
   std::vector<double> ret;
   ret.push_back(type);
   ret.push_back(degrees);
   return ret;
}

std::string FluxMgr::writeXmlFile(const std::vector<std::string>& fileList) {
/** purpose: creates a document of the form

  <?xml version='1.0' ?>
  <!DOCTYPE source_library SYSTEM "d:\users\burnett\pdr_v7r1c\flux\v5r3/xml/source.dtd" [
  <!ENTITY librarya SYSTEM "d:\users\burnett\pdr_v7r1c\flux\v5r3/xml/user_library.xml" >
  <!ENTITY libraryb SYSTEM "d:\users\burnett\pdr_v7r1c\flux\v5r3/xml/source_library.xml" >
  ]>
  <source_library>
  &librarya;
  &libraryb;
  </source_library>
  
    */
    std::stringstream fileString;
    // Unique tag to add to ENTITY elements in the DTD.
    char libchar = 'a';
    std::string inFileName;
    
    std::vector<std::string>::const_iterator iter = fileList.begin();
    
    //the default DTD file
    inFileName=m_dtd;
    //replace $(FLUXROOT) by its system variable
    xml::IFile::extractEnvVar(&inFileName);
    
    //this stuff goes in the beginnning of the XML file to be read into the parser
    fileString << "<?xml version='1.0' ?>" << std::endl << "<!DOCTYPE source_library" 
        << " SYSTEM " << '"' << inFileName << '"' << " [" << std::endl;
    
    //as long as there are files in the file list...
    for (;iter != fileList.end(); iter++) {
        
        // get the file name, and evaluate any system variables in it
        inFileName=(*iter).c_str();
        xml::IFile::extractEnvVar(&inFileName);
        
        //then make an ENTITY entry specifying where the file is
        fileString << "<!ENTITY " << "library" << libchar << " SYSTEM " << '"' 
            << inFileName << "\" >" << std::endl;      
        libchar++;
    }
    
    fileString << "]>" << std::endl << "<source_library>" << std::endl;
    iter = fileList.begin();
    libchar = 'a';
    
    //as long as there are files in the file list...
    for (;iter != fileList.end(); iter++) {
        // add a reference to the file name
        fileString << "&library" << libchar << ";" << std::endl;       
        libchar++;
    }
    
    fileString << "</source_library>" << '\0';
    return fileString.str();
    
}
