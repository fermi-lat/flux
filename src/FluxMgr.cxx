/** @file FluxMgr.cpp
@brief Implementation of FluxMgr
*/

#include "flux/FluxMgr.h"
#include "flux/EventSource.h"
#include "flux/SpectrumFactoryTable.h"
#include "astro/GPS.h"
#include "flux/FluxException.h" // defines FATAL_MACRO
#include "flux/CompositeSource.h"

#include "facilities/Util.h"     // for expandEnvVar
#include "facilities/commonUtilities.h"

#include "astro/PointingTransform.h"

// RapidXML framework
#include "xmlBase/rapidxml.hpp"
#include "xmlBase/rapidxml_error_framework.hpp"
#include "xmlBase/safe_xml_parser.hpp"

#include <sstream>
#include <map>
#include <iostream>
#include <iterator>
#include <stdexcept>
#include <string>
#include <fstream>
#include <cstring>
#include <algorithm> 

#define DECLARE_SPECTRUM(x)   extern const ISpectrumFactory& x##Factory; x##Factory.addRef();

using astro::GPS;
using namespace xml_framework;

// Helper function implementations
std::string FluxMgr::getAttribute(XmlNode* node, const char* name) {
    if (!node) return "";
    auto* attr = node->first_attribute(name);
    return attr ? std::string(attr->value()) : "";
}

double FluxMgr::getDoubleAttribute(XmlNode* node, const char* name) {
    std::string val = getAttribute(node, name);
    if (val.empty()) return 0.0;
    try {
        return std::stod(val);
    } catch (...) {
        return 0.0;
    }
}

XmlNode* FluxMgr::findFirstChildByName(XmlNode* parent, const char* name) {
    if (!parent) return nullptr;
    // Handle wildcard "*" to get first child element
    if (name && std::strcmp(name, "*") == 0) {
        return getFirstChildElement(parent);
    }
    return parent->first_node(name);
}

XmlNode* FluxMgr::getFirstChildElement(XmlNode* parent) {
    if (!parent) return nullptr;
    for (auto* child = parent->first_node(); child; child = child->next_sibling()) {
        if (child->type() == rapidxml::node_element) {
            return child;
        }
    }
    return nullptr;
}

XmlNode* FluxMgr::getSiblingElement(XmlNode* node) {
    if (!node) return nullptr;
    for (auto* sibling = node->next_sibling(); sibling; sibling = sibling->next_sibling()) {
        if (sibling->type() == rapidxml::node_element) {
            return sibling;
        }
    }
    return nullptr;
}

std::string FluxMgr::getTagName(XmlNode* node) {
    if (!node || !node->name()) return "";
    return std::string(node->name(), node->name_size());
}

FluxMgr::FluxMgr(const std::vector<std::string>& fileList, std::string dtdname)
: m_dtd(dtdname.empty() ? 
        facilities::commonUtilities::joinPath(
            facilities::commonUtilities::getXmlPath("flux"), "source.dtd") 
        : dtdname)
{
    if (fileList.empty()) {
        defaultFile();
    } else {		
        init(fileList);		
    }	
}

void FluxMgr::defaultFile() {
    // Purpose: to set the default xml file and initialize the package to use it.
    std::vector<std::string> input;

    // must find the source_library.xml file.
    input.push_back(facilities::commonUtilities::joinPath(
        facilities::commonUtilities::getXmlPath("flux"), "source_library.xml"));
    init(input);
}

void FluxMgr::init(const std::vector<std::string>& fileList) {
    // Parse each XML file individually and collect sources
    for (const auto& filename : fileList) {
        std::string expandedFilename = filename;
        facilities::Util::expandEnvVar(&expandedFilename);
        
        // Read file into buffer
        std::ifstream file(expandedFilename, std::ios::binary | std::ios::ate);
        if (!file.is_open()) {
            std::cerr << "Warning: Could not open file: " << expandedFilename << std::endl;
            continue;
        }
        
        auto size = file.tellg();
        file.seekg(0, std::ios::beg);
        
        // Create buffer for this file (RapidXML modifies the buffer during parsing)
        m_xmlBuffers.emplace_back(static_cast<size_t>(size) + 1);
        auto& buffer = m_xmlBuffers.back();
        
        if (!file.read(buffer.data(), size)) {
            std::cerr << "Warning: Failed to read file: " << expandedFilename << std::endl;
            m_xmlBuffers.pop_back();
            continue;
        }
        buffer[static_cast<size_t>(size)] = '\0';
        
        // Create and parse document
        auto doc = std::make_unique<XmlDocument>();
        try {
            doc->parse<rapidxml::parse_default>(buffer.data());
        } catch (const rapidxml::parse_error& e) {
            std::cerr << "XML parse error in " << expandedFilename << ": " << e.what() << std::endl;
            m_xmlBuffers.pop_back();
            continue;
        }
        
        // Iterate through ALL source_library elements in the document
        // (handles XML files with multiple <source_library> sections)
        for (XmlNode* root = doc->first_node("source_library"); 
             root; 
             root = root->next_sibling("source_library")) {
            
            // Set s_library to the first valid root if not already set
            if (!s_library) {
                s_library = root;
            }
            
            // Get title attribute for this library
            std::string libraryTitle = getAttribute(root, "title");
            if (libraryTitle.empty()) {
                libraryTitle = expandedFilename;
            }
            
            // Iterate through sources in this library
            for (XmlNode* source = root->first_node("source"); 
                 source; 
                 source = source->next_sibling("source")) {
                std::string sourceName = getAttribute(source, "name");
                if (!sourceName.empty()) {
                    m_sources[sourceName] = std::make_pair(source, libraryTitle);
                }
            }
        }
        
        // Fallback: if no source_library found, try the document root directly
        // (for backward compatibility with simpler XML structures)
        if (m_sources.empty()) {
            XmlNode* root = doc->first_node();
            if (root) {
                if (!s_library) {
                    s_library = root;
                }
                
                std::string libraryTitle = getAttribute(root, "title");
                if (libraryTitle.empty()) {
                    libraryTitle = expandedFilename;
                }
                
                for (XmlNode* source = root->first_node("source"); 
                     source; 
                     source = source->next_sibling("source")) {
                    std::string sourceName = getAttribute(source, "name");
                    if (!sourceName.empty()) {
                        m_sources[sourceName] = std::make_pair(source, libraryTitle);
                    }
                }
            }
        }
        
        m_documents.push_back(std::move(doc));
    }
    
    // Register locally defined spectra
    DECLARE_SPECTRUM(TimeCandle);
    DECLARE_SPECTRUM(FileSource);
    DECLARE_SPECTRUM(SurfaceMuons);
    DECLARE_SPECTRUM(VdgGamma);
    DECLARE_SPECTRUM(Earth);
}

FluxMgr::~FluxMgr() {
    // unique_ptr handles cleanup automatically
}

EventSource* FluxMgr::source(std::string name) {
    // Purpose: to return a pointer to a source, referenced by name.
    auto it = m_sources.find(name);
    if (it == m_sources.end()) {
        return nullptr;
    }
    return getSourceFromXML(it->second.first);
}

EventSource* FluxMgr::compositeSource(std::vector<std::string> names) {
    // Purpose: to return a pointer to a source, referenced by a list of names.
    CompositeSource* comp = new CompositeSource();
    
    for (const auto& name : names) {
        if (m_sources.find(name) == m_sources.end()) {
            std::cerr << "Unrecognized source: " << name << std::endl;
            std::cerr << "Known names: " << std::endl;
            std::list<std::string> known(sourceList());
            std::copy(known.begin(), known.end(), 
                      std::ostream_iterator<std::string>(std::cerr, ", "));
            FATAL_MACRO("Unrecognized source name " + name);
            delete comp;
            return nullptr;
        }
        comp->addSource(getSourceFromXML(m_sources[name].first));
    }
    return comp;
}

EventSource* FluxMgr::getSourceFromXML(XmlNode* src) {
    // Purpose: sourceFromXML - create a new EventSource from a DOM element
    if (!src) {
        FATAL_MACRO("Null XML source element");
        return nullptr;
    }
    
    XmlNode* childNode = getFirstChildElement(src);
    if (childNode == nullptr) {
        // no child node: expect to find the name defined.
        return new FluxSource(src);
    }
    
    // Check if this is a nested source reference
    std::string childTag = getTagName(childNode);
    if (childTag == "nestedSource") {
        std::string nestedName = getAttribute(childNode, "sourceRef");
        auto it = m_sources.find(nestedName);
        if (it == m_sources.end()) {
            FATAL_MACRO("Nested source not found: " + nestedName);
            return nullptr;
        }
        return getSourceFromXML(it->second.first);
    }
    
    return new FluxSource(src);
}

XmlNode* FluxMgr::getLibrarySource(const std::string& id) {
    // Purpose: source library lookup by ID
    auto it = m_sources.find(id);
    return (it != m_sources.end()) ? it->second.first : nullptr;
}

std::list<std::string> FluxMgr::sourceList() const {
    std::list<std::string> s;
    for (const auto& [name, pair] : m_sources) {
        s.push_back(name);
    }
    return s;
}

std::vector<std::pair<std::string, std::list<std::string>>> FluxMgr::sourceOriginList() const {
    std::vector<std::pair<std::string, std::list<std::string>>> originList;
    
    for (const auto& [name, pair] : m_sources) {
        const std::string& filename = pair.second;
        
        // Find if we already have an entry for this file
        auto it = std::find_if(originList.begin(), originList.end(),
            [&filename](const auto& entry) { return entry.first == filename; });
        
        if (it != originList.end()) {
            it->second.push_back(name);
        } else {
            originList.emplace_back(filename, std::list<std::string>{name});
        }
    }
    return originList;
}

double FluxMgr::time() const {
    return GPS::instance()->time();
}

void FluxMgr::setTime(double newtime) {
    GPS::instance()->time(newtime);
}

void FluxMgr::synch() {
    GPS::instance()->synch();
}

void FluxMgr::pass(double t) {
    GPS::instance()->pass(t);
}

std::pair<double, double> FluxMgr::location() {
    return std::make_pair(GPS::instance()->lat(), GPS::instance()->lon());
}

CLHEP::HepRotation FluxMgr::transformToGlast(double seconds, GPS::CoordSystem index) {
    return GPS::instance()->transformToGlast(seconds, index);
}

std::vector<double> FluxMgr::setRockType(GPS::RockType rockType, double rockAngle) {
    int type = GPS::instance()->setRockType(rockType);
    double degrees = GPS::instance()->rockingDegrees(rockAngle);
    return {static_cast<double>(type), degrees};
}

void FluxMgr::setAlignmentRotation(double qx, double qy, double qz, bool misalign) {
    CLHEP::HepRotation R(
        CLHEP::HepRotationX(qx * M_PI / 180) * 
        CLHEP::HepRotationY(qy * M_PI / 180) * 
        CLHEP::HepRotationZ(qz * M_PI / 180));
    
    if (misalign) {
        EventSource::setAlignmentRotation(R);
    } else {
        GPS::instance()->setAlignmentRotation(R);
    }
}

int FluxMgr::setIdOffset(int id) {
    int last = EventSource::s_id_offset;
    EventSource::s_id_offset = id;
    return last;
}

void FluxMgr::setFilterCone(double ra, double dec, double radius) {
    EventSource::s_cone = {ra, dec, radius};
}

void FluxMgr::setArea(double area) {
   EventSource::totalArea(area);
}

void FluxMgr::setPointingHistoryFile(std::string fileName) {
    GPS::instance()->setPointingHistoryFile(fileName);
}

void FluxMgr::addFactory(std::string name, const ISpectrumFactory* factory) {
    SpectrumFactoryTable::instance()->addFactory(name, factory);
}

void FluxMgr::setExpansion(double p) {
    GPS::instance()->expansion(p);
}

/// generate some test output
void FluxMgr::test(std::ostream& cout, std::string source_name, int count)
{   
    using astro::GPS;
    EventSource* e = source(source_name);
    if (e == nullptr) {
        throw std::invalid_argument(std::string("Did not find source ") + source_name);
    }
    setExpansion(1.);
    double time = 0.;

    const int howMany = e->howManySources();
    std::map<int, int> counts;

    cout << "running source: " << e->fullTitle() << std::endl;
    cout << " Total rate is: " << e->rate(time) << " Hz into " << e->totalArea() << " m^2" << std::endl;
    cout << " there are " << howMany << " Sources total..." << std::endl;
    cout << "    Generating " << count << " trials " << std::endl;
    cout << " --------------------------------" << std::endl;

    EventSource* f = nullptr;
    double totalinterval = 0.;
    
    for (int i = 0; i < count; ++i) {
        f = e->event(time);

        if (!f->enabled()) {
            std::cout << "Source turned off at time " << time << std::endl;
            break;
        }
        
        double interval = e->interval();

        // Increment the "elapsed" time and the "orbital" time,
        // just as is done in flux. NOTE: this is important for the operation 
        // of fluxsource, and is expected.
        time += interval;
        pass(interval);
        
        int sourceNumber = e->numSource();
        if (sourceNumber == -1) {
            counts[0]++;
        } else {
            counts[sourceNumber]++;
        }

        totalinterval += interval;
        
        cout << f->particleName()
             << "(" << f->energy() << " MeV)"
             << ", Launch: " << f->launchPoint() 
             << ", Dir " << f->launchDir() 
             << ", Flux=" << f->flux(time) 
             << ", Interval=" << interval;
        
        if (sourceNumber != -1) {
            cout << ", SourceID: " << sourceNumber;
        }
        
        cout << "\tElapsed time= " << totalinterval 
             << std::endl;
    }
    
    cout << "------------------------------------------------------" << std::endl;

    cout << std::endl << "Average Interval=" << totalinterval / count << " , "
         << "Average rate = " << count / totalinterval << std::endl;

    cout << "Source Statistics: " << std::endl;
    for (const auto& [sourceId, eventCount] : counts) {
        cout << "source #" << sourceId << ": " << eventCount << " events counted." << std::endl;
    }
}
