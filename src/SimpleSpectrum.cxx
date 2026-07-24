/** @file SimpleSpectrum.cpp
    @brief definition of SimpleSpectrum
*/

#include "flux/SimpleSpectrum.h"
#include "flux/SpectrumFactory.h"
#include "flux/FluxException.h" // for FATAL_MACRO

// RapidXML includes
#include "xmlBase/rapidxml.hpp"

#include "facilities/Util.h"

#include <cstdlib>
#include <utility>
#include <iostream>
#include <sstream>
#include <cmath>
#include <map>
#include <stdexcept>
#include <cstring>

static SpectrumFactory<SimpleSpectrum> factory;

namespace {
    // useful utility functions

    // differential rate: return energy distributed as e**-gamma between e1 and e2, 
    // if r is uniform from 0 to 1
    double power_law(double r, double e1, double e2, double gamma) {
        double e = gamma == 1
            ? e1 * exp(r * log(e2 / e1))
            : e1 * exp(log(1.0 - r * (1. - pow(e2 / e1, 1 - gamma))) / (1 - gamma));
        return e;
    }
    
    // integral of e**-gamma from e1 to e2
    double total_rate(double e1, double e2, double gamma) {
        return gamma == 1
            ? log(e2 / e1)
            : (pow(e1, 1 - gamma) - pow(e2, 1 - gamma)) / (gamma - 1);
    }
    
    ///@class ParMap 
    ///@brief local analysis of keyword string
    class ParMap {
    public:
        ParMap(std::string paramString) {
            facilities::Util::keyValueTokenize(paramString, ",", m_tokenMap);
        }
        
        double value(std::string name) const {
            auto it = m_tokenMap.find(name);
            if (it == m_tokenMap.end()) {
                throw std::invalid_argument("SimpleSpectrum: keyword " + name + " not found");
            }
            return std::atof(it->second.c_str());
        }
        
    private:
        std::map<std::string, std::string> m_tokenMap;
    };
}

// Helper function implementations
std::string SimpleSpectrum::getAttribute(XmlNode* node, const char* name) {
    if (!node) return "";
    auto* attr = node->first_attribute(name);
    return attr ? std::string(attr->value()) : "";
}

double SimpleSpectrum::getDoubleAttribute(XmlNode* node, const char* name) {
    std::string val = getAttribute(node, name);
    if (val.empty()) return 0.0;
    try {
        return std::stod(val);
    } catch (...) {
        return 0.0;
    }
}

XmlNode* SimpleSpectrum::findFirstChildByName(XmlNode* parent, const char* name) {
    if (!parent) return nullptr;
    // Handle wildcard "*" to get first child element
    if (name && std::strcmp(name, "*") == 0) {
        for (auto* child = parent->first_node(); child; child = child->next_sibling()) {
            if (child->type() == rapidxml::node_element) {
                return child;
            }
        }
        return nullptr;
    }
    return parent->first_node(name);
}

std::string SimpleSpectrum::getTagName(XmlNode* node) {
    if (!node || !node->name()) return "";
    return std::string(node->name(), node->name_size());
}

// implement simple broken power law
SimpleSpectrum::SimpleSpectrum(const std::string& paramString)
    : m_name("gamma")
    , m_E0(10)
    , m_index(2.0)
    , m_index2(2.0)
    , m_ebreak(0)
    , m_emax(200000)
    , m_useGeV(false)
    , m_a(1.0)
{
    ParMap parmap(paramString);
    
    try {
        m_E0 = parmap.value("emin");
    } catch (...) { }
    
    try {
        m_emax = parmap.value("emax");
    } catch (...) {}
    
    try {
        m_index = parmap.value("gamma");
    } catch (...) {}
    
    try {
        m_index2 = parmap.value("gamma2");
        m_ebreak = parmap.value("ebreak");
    } catch (...) {}

    setup_power_law();
}

SimpleSpectrum::SimpleSpectrum(XmlNode* xelem, bool useGeV)
    : m_useGeV(useGeV)
    , m_E0(10)
    , m_index(2.0)
    , m_index2(2.0)
    , m_ebreak(0)
    , m_emax(200000)
    , m_a(1.0)
{
    m_name = getAttribute(xelem, "name");
    
    XmlNode* spectrum = findFirstChildByName(xelem, "*");
    
    std::string tagName = getTagName(spectrum);
    
    if (tagName == "power_law") {
        m_E0 = getDoubleAttribute(spectrum, "emin");
        m_emax = getDoubleAttribute(spectrum, "emax");
        m_index = getDoubleAttribute(spectrum, "gamma");
        m_ebreak = getDoubleAttribute(spectrum, "ebreak");
        m_index2 = getDoubleAttribute(spectrum, "gamma2");
        setup_power_law();
    }
    else if (tagName == "energy") {
        // single energy: no interpolation
        m_emax = m_E0 = getDoubleAttribute(spectrum, "e");
    }
    else if (tagName == "exponential") {
        m_E0 = getDoubleAttribute(spectrum, "exponential");
        m_index = getDoubleAttribute(spectrum, "exponent");
        m_emax = 100.0;
        m_index = 0.0;
        FATAL_MACRO("exponential spectral component not implemented yet");
    }
    else {
        std::cerr << "Unknown name: " << m_name << std::endl;
        FATAL_MACRO("Unknown particle spectrum!");
    }
}

void SimpleSpectrum::setup_power_law() {
    if (m_ebreak == 0) {
        m_ebreak = m_emax;
        m_a = 1.0; 
    } else {
        // calculate relative part of spectrum for lower
        double a1 = total_rate(m_E0, m_ebreak, m_index);
        double a2 = pow(m_ebreak, m_index2 - m_index) * total_rate(m_ebreak, m_emax, m_index2);
        m_a = a1 / (a1 + a2);
    }
}

std::string SimpleSpectrum::title() const {
    std::stringstream s;
    s << particleName() << '(' << m_E0 << (m_useGeV ? " GeV" : " MeV");
    if (m_index >= 1) s << ',' << m_index;
    if (m_ebreak != 0) s << "," << m_ebreak << "," << m_index2;
    s << ")";
    return s.str();
}

float SimpleSpectrum::operator()(float f) {
    if (m_emax == m_E0) return m_E0;
    
    float energy;
    if (f < m_a) {
        // single power law, or first segment
        energy = power_law(f / m_a, m_E0, m_ebreak, m_index);
    } else {
        // break in the power law above the break
        energy = power_law((f - m_a) / (1 - m_a), m_ebreak, m_emax, m_index2);
    }
    return energy;
}

const char* SimpleSpectrum::particleName() const {
    return m_name.c_str();
}

float SimpleSpectrum::parseParamList(std::string input, int index) {
    std::vector<float> output;
    int i = 0;
    for (; !input.empty() && i != std::string::npos;) {
        float f = std::atof(input.c_str());
        output.push_back(f);
        i = input.find_first_of(",");
        input = input.substr(i + 1);
    }
    
    if (static_cast<size_t>(index) < output.size()) {
        return output.at(index);
    } else {
        return 0;
    }
}
