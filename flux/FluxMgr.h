/** @file FluxMgr.h
    @brief declaration of FluxMgr

 $Header: /nfs/slac/g/glast/ground/cvs/FluxSvc/src/FluxMgr.h,v 1.20 2003/03/21 19:14:37 jrb Exp $

  */
#ifndef FLUX_MGR_H
#define FLUX_MGR_H

/** 
* \class FluxMgr
*
* \brief The point of entry for interfacing with the flux package.
* holds methods for creating sources and sending new particles, 
* and methods for interfacing with the satellite position, and 
* setting the position variables. It is instantiated with
* the names of the xml files to be used as input to the xml parser.
* 
* $Header $
*/

#include "GPS.h"

#include "FluxSource.h"

#include <xercesc/dom/DOM_Document.hpp>
#include <xercesc/dom/DOM_Element.hpp>
#include "xml/XmlParser.h"
#include "ISpectrumFactory.h"
#include <map>
#include <list>
#include <string>

class FluxMgr 
{
    
public:
    
    /// ctor for multiple XML documents
    FluxMgr(const std::vector<std::string>& fileList, std::string dtd="");
    
    ~FluxMgr();
    
    /// create and return a source by name.
    EventSource* source(std::string name);
    
    /// access to the source list
    std::list<std::string> sourceList() const;
    
    /// set the target area
    void setArea(double area);
    
    /// generate some test output
    void test(std::ostream& out, std::string source_name, int count);
    
    /// set the angular (off-zenith) values of the GLAST satellite
    void setExplicitRockingAngles(std::pair<double,double> ang);

    /// get the angular values of the satellite
    std::pair<double,double> getExplicitRockingAngles();

    ///this should return the source file names, along with the contained sources.
    std::vector<std::pair< std::string ,std::list<std::string> > > sourceOriginList() const;
    
    void addFactory(std::string name, const ISpectrumFactory* factory );
    
    /// set the expansion factor for the orbit (-1) = random
    void setExpansion (double p);
    
    /// pass a specific amount of time
    void pass ( double t);
    
    /// Get the time as held by GPS
    GPStime time () const;
    
    /// synch satellite location with current time
    void synch ();
    
    /// set the sample interval
    void sampleintvl ( /*GPStime*/double t );
    
    /// get the current satellite location
    std::pair<double,double> location();
    
    
    ///get the transformation matrix due to orientation of the Galaxy
    HepRotation CELTransform(double time);
    
    ///get the transformation matrix due to orientation of the spacecraft.
    HepRotation orientTransform(double time);
    
    ///this transforms glast-local (cartesian) vectors into galactic (cartesian) vectors
    HepRotation FluxMgr::transformGlastToGalactic(double time);

    ///this sets the rocking mode in GPS.
    std::vector<double> setRockType(GPS::RockType rockType, double rockAngle);
    std::vector<double> setRockType(int rockType, double rockAngle);
private:
    
    /// source library lookup.  Each source is uniquely identified
    /// by its "name" attribute because "name" is of type ID
    DOM_Element  getLibrarySource(const DOMString& id);
    
    
    
    void defaultFile();
    void init(const std::vector<std::string>& fileList);
    
    EventSource* getSourceFromXML(const DOM_Element& src);
    
    DOM_Document m_library_doc;
    
    DOM_Element     s_library;
    
    std::vector<DOM_Document> m_library_doclist;
    
    std::vector<DOM_Element>     s_librarylist;
    
    /// list of sources for easy lookup
    //std::map<std::string, DOM_Element > m_sources;
    std::map<std::string, std::pair<DOM_Element,std::string> > m_sources;

    /// internal routine that creates the document
    std::string  writeXmlFile( const std::vector<std::string>& fileList);
    
    /// filename for dtd
    std::string m_dtd;
};
#endif
