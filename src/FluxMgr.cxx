/** @file FluxMgr.cxx
@brief Implementation of FluxMgr

$Header: /nfs/slac/g/glast/ground/cvs/flux/src/FluxMgr.cxx,v 1.18 2004/09/29 06:52:28 cohen Exp $
*/

#include "flux/FluxMgr.h"
#include "flux/EventSource.h"
#include "flux/SpectrumFactoryTable.h"
#include "astro/GPS.h"
#include "flux/FluxException.h" // defines FATAL_MACRO
#include "flux/CompositeSource.h"

// #include <xercesc/dom/DOMDocument.hpp>   already included by .h file
// #include <xercesc/dom/DOMElement.hpp>
#include "xml/Dom.h"
#include "facilities/Util.h"     // for expandEnvVar

#include "astro/PointingTransform.h"

#include <sstream>
#include <map>
#include <stdexcept>
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

    if (m_library_doc == 0) {
        FATAL_MACRO("Parse error: processing the document" << std::endl
            << xmlFileIn << std::endl);
        return;
    }

    // Root element is of type source_library.  Content is
    // one or more source elements.

    s_library = m_library_doc->getDocumentElement();

    // loop through the source elements to create a map of names, DOMElements
    if (s_library != 0) {

        DOMElement* child = xml::Dom::getFirstChildElement(s_library);
        DOMElement* toplevel = xml::Dom::getFirstChildElement(s_library);

        while (child != 0) {
            while (!(xml::Dom::hasAttribute(child, "name"))  )
            {
                s_library = child;
                child = xml::Dom::getFirstChildElement(s_library);
            }

            while (child != 0) {
                std::string name = xml::Dom::getAttribute(child, "name");
                std::string parentfilename = xml::Dom::getAttribute(toplevel, "title");
                m_sources[name]=std::make_pair<DOMElement*,std::string>(child,parentfilename);
                child = xml::Dom::getSiblingElement(child);
            }

            child = xml::Dom::getSiblingElement(toplevel);
            toplevel=child;
        }

    }
    // these are the locally defined spectra that we want to make available
    DECLARE_SPECTRUM( FILESpectrum);
    DECLARE_SPECTRUM( TimeCandle);

    DECLARE_SPECTRUM( SurfaceMuons);
    DECLARE_SPECTRUM( GalPulsars);
    DECLARE_SPECTRUM( VdgGamma);

    // these are deprecated, will be replaced by Hiroshima group
    DECLARE_SPECTRUM( AlbedoPSpectrum);
    DECLARE_SPECTRUM( CHIMESpectrum );
    DECLARE_SPECTRUM( GalElSpectrum);
    DECLARE_SPECTRUM( MapSpectrum);
    DECLARE_SPECTRUM( AGNSpectrum);

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

EventSource* FluxMgr::compositeSource(std::vector<std::string> names)
{
    //Purpose: to return a pointer to a source, referenced by a list of names.
    //Input: the names of the desired sources.

    CompositeSource* comp = new CompositeSource();
    for( std::vector<std::string>::const_iterator it= names.begin(); it!=names.end(); ++it){
        const std::string& name = *it;
        if( m_sources.find(name)==m_sources.end() ) {
            delete comp;
            return 0;
        }
        comp->addSource(getSourceFromXML(m_sources[name].first));
    }
    return comp;
}

EventSource*  FluxMgr::getSourceFromXML(const DOMElement* src)
{
    //Purpose: sourceFromXML - create a new EventSource from a DOM element
    //instantiated, e.g., from a description in source_library.xml
    //Input:  the element holding particle information.

    using XERCES_CPP_NAMESPACE_QUALIFIER DOMNode;

    DOMNode*    childNode = src->getFirstChild();
    if (childNode == 0) {
        /*
        FATAL_MACRO("Improperly formed XML event source");
        return 0;
        */
        // no child node: expect to find the name defined.
        return  new FluxSource(src);
    }

    DOMElement* sname = xml::Dom::getFirstChildElement(src);
    if (sname == 0 ) {
        FATAL_MACRO("Improperly formed XML event source");
        return 0;
    }
    // If we got here, should have legit child element
    if (xml::Dom::checkTagName(sname, "spectrum") )
    {
        return  new FluxSource(src);
    }
    else if (xml::Dom::checkTagName(sname, "nestedSource"))
    {

        // Search for and process immediate child elements.  All must
        // be of type "nestedSource".  There may be more than one.
        // Content model for nestedSource is EMPTY, so can omit check
        // for that in the code

        CompositeSource* cs;
        cs = new CompositeSource();
        do { 
            DOMElement* selem = 
                getLibrarySource(xml::Dom::getAttribute(sname, "sourceRef"));

            if (selem == 0) {
                FATAL_MACRO("source name" << 
                    xml::Dom::getAttribute(sname, "sourceRef") <<
                    "' not in source library");
            }
            cs->addSource(getSourceFromXML(selem)); 
            sname = xml::Dom::getSiblingElement(sname);
        } 
        while (sname != 0 );
        return cs;
    }
    else {
        FATAL_MACRO("Unexpected element: " << 
            xml::Dom::getTagName(sname) );
    }
    return 0;
}



DOMElement*    FluxMgr::getLibrarySource(const std::string& id)
{
    //Purpose: source library lookup.  Each source is uniquely identified
    // by its "name" attribute because "name" is of type ID

    // quit if the library was unitialized
    if (s_library == 0 ) return 0; 

    return xml::Dom::getElementById(m_library_doc, id);
}

std::list<std::string> FluxMgr::sourceList() const
{
    std::list<std::string> s;
    for( std::map<std::string, std::pair<DOMElement*,std::string> >::const_iterator it = m_sources.begin();
        it != m_sources.end();
        ++it){
            s.push_back((*it).first);
        }
        return s;
}

std::vector<std::pair< std::string ,std::list<std::string> > > FluxMgr::sourceOriginList() const
{
    std::vector<std::pair< std::string ,std::list<std::string> > > originList;
    for( std::map<std::string, std::pair<DOMElement*,std::string> >::const_iterator it = m_sources.begin();
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
    if (e==0) {
        throw std::invalid_argument(std::string("Did not find source ")+source_name);
    }
    setExpansion(1.);
    double time=0.;

    const int howMany = e->howManySources();
    std::map<int,int> counts;


    cout << "running source: " << e->fullTitle() << std::endl;
    cout << " Total rate is: " << e->rate(time) << " Hz into " << e->totalArea() << " m^2" << std::endl;
    //cout << "LaunchType" << f->retLaunch() << "Pointtype" << f->retPoint() <<std::endl;
    cout << " there are " << howMany << " Sources total..." << std::endl;
    cout << "    Generating " << count << " trials " << std::endl;
    cout << " --------------------------------" << std::endl;

    //testing rotateangles function
    GPS::instance()->rotateAngles(std::make_pair<double,double>(0.0,0.0));
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
            << "(" << f->energy()<< " MeV)"
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
    for( std::map<int,int>::const_iterator q=counts.begin() ; q!= counts.end() ; ++q){
        cout << "source #" << q->first << ": " << q->second << " events counted." << std::endl;
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

//get the transformtation matrix - the rest of these functions are now deprecated
HepRotation FluxMgr::transformToGlast(double seconds,GPS::CoordSystem index){
    return GPS::instance()->transformToGlast(seconds, index);
}

//get the transformation matrix.
HepRotation FluxMgr::CELTransform(double time){
    return GPS::instance()->CELTransform(time);
}

//get the transformation matrix.
HepRotation FluxMgr::orientTransform(double time){
    //make the transformtion that turns zenith coordinates into local coordinates.
    HepRotation ret;
    ret = GPS::instance()->transformToGlast(time,GPS::ZENITH);
    return ret;
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
    facilities::Util::expandEnvVar(&inFileName);

    //this stuff goes in the beginnning of the XML file to be read into the parser
    fileString << "<?xml version='1.0' ?>" << std::endl << "<!DOCTYPE source_library" 
        << " SYSTEM " << '"' << inFileName << '"' << " [" << std::endl;

    //as long as there are files in the file list...
    for (;iter != fileList.end(); iter++) {

        // get the file name, and evaluate any system variables in it
        inFileName=(*iter).c_str();
        facilities::Util::expandEnvVar(&inFileName);

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
