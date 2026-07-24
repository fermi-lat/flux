/** @file FluxMgr.h
    @brief declaration of FluxMgr
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
*/

#include "astro/GPS.h"
#include "FluxSource.h"

// RapidXML-based framework includes
#include "xmlBase/rapidxml.hpp"
#include "xmlBase/rapidxml_error_framework.hpp"
#include "xmlBase/safe_xml_parser.hpp"
#include "xmlBase/xml_result.hpp"

#include "ISpectrumFactory.h"
#include <map>
#include <list>
#include <string>
#include <vector>
#include <memory>

// Type aliases for RapidXML types
using XmlNode = rapidxml::xml_node<char>;
using XmlDocument = rapidxml::xml_document<char>;

class FluxMgr 
{
    
public:
    
    /// ctor for multiple XML documents
    FluxMgr(const std::vector<std::string>& fileList, std::string dtd="");
    
    ~FluxMgr();
    
    /// create and return a source by name.
    EventSource* source(std::string name);

    /// create a composite source from the list of names
    EventSource* compositeSource(std::vector<std::string> names);

    /// access to the source list
    std::list<std::string> sourceList() const;
    
    /// set the target area
    void setArea(double area);
    
    /// generate some test output
    void test(std::ostream& out, std::string source_name, int count);

    /// set the desired pointing history file to use:
    void setPointingHistoryFile(std::string fileName);

    ///this should return the source file names, along with the contained sources.
    std::vector<std::pair< std::string ,std::list<std::string> > > sourceOriginList() const;
    
    void addFactory(std::string name, const ISpectrumFactory* factory );
    
    /// set the expansion factor for the orbit (-1) = random
    void setExpansion (double p);

    /// pass a specific amount of time
    void pass ( double t);

    /// Get the time as held by GPS
    double time () const;

    /// Set the time
    void setTime(double newtime);

    /// synch satellite location with current time
    void synch ();
    
    /// set the sample interval
    void sampleintvl ( /*GPStime*/double t );
    
    /// get the current satellite location
    std::pair<double,double> location();
    
    CLHEP::HepRotation transformToGlast(double seconds, astro::GPS::CoordSystem index);

    ///this sets the rocking mode in GPS.
    std::vector<double> setRockType(astro::GPS::RockType rockType, double rockAngle);
    
    /// Set an alignment rotation:
    /// @param qx,qy, qz rotation angles (degrees) about coordinate axes -- assume small, < 1 deg
    /// @param misalign if true, apply to incoming; if false, apply as correction to coordinate transformation
    void setAlignmentRotation(double qx, double qy, double qz, bool misalign);

    /// set an offset for generating source id numbers, return previous value
    int setIdOffset(int id);

    /** set a cone to filter incoming (galactic) data
    @param ra,dec center of cone, equatorial coords in degrees
    @param radius radius of cone, degrees
    */
    void setFilterCone(double ra, double dec, double radius);

private:
    
    /// source library lookup.  Each source is uniquely identified
    /// by its "name" attribute because "name" is of type ID
    XmlNode* getLibrarySource(const std::string& id);
    
    void defaultFile();
    void init(const std::vector<std::string>& fileList);
    
    EventSource* getSourceFromXML(XmlNode* src);
    
    /// Owned document storage - each file gets its own document
    std::vector<std::unique_ptr<XmlDocument>> m_documents;
    
    /// Buffer storage for parsed XML content (RapidXML requires persistent buffers)
    std::vector<std::vector<char>> m_xmlBuffers;
    
    /// Root library element pointer
    XmlNode* s_library = nullptr;
    
    /// list of sources for easy lookup
    std::map<std::string, std::pair<XmlNode*, std::string>> m_sources;

    /// filename for dtd (kept for compatibility, but DTD validation not used with RapidXML)
    std::string m_dtd;
    
    // Helper functions for XML operations
    static std::string getAttribute(XmlNode* node, const char* name);
    static double getDoubleAttribute(XmlNode* node, const char* name);
    static XmlNode* findFirstChildByName(XmlNode* parent, const char* name);
    static XmlNode* getFirstChildElement(XmlNode* parent);
    static XmlNode* getSiblingElement(XmlNode* node);
    static std::string getTagName(XmlNode* node);
};
#endif
