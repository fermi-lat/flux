/**
 * @file FileSource.cxx
 * @brief Read in the incident particle properties from a file.
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/flux/src/FileSource.cxx,v 1.9 2006/11/05 20:08:42 burnett Exp $
 */

#include <cstdlib>
#include <fstream>
#include <map>

#include "facilities/Util.h"

#include "flux/SpectrumFactory.h"
#include "flux/EventSource.h"

#include "FileSource.h"

static SpectrumFactory<FileSource> factory;
const ISpectrumFactory& FileSourceFactory = factory;

namespace {
   void readLines(std::string inputFile, 
                  std::vector<std::string> &lines,
                  const std::string &skip) {
      facilities::Util::expandEnvVar(&inputFile);
      std::ifstream file(inputFile.c_str());
      lines.clear();
      std::string line;
      while (std::getline(file, line, '\n')) {
         if (line != "" && line != " "             //skip (most) blank lines 
             && line.find_first_of(skip) != 0) {   //and commented lines
            lines.push_back(line);
         }
      }
   }
}

FileSource::FileSource(const std::string & params) 
   : m_launchDirection(0), m_launchPoint(0), m_currentLine(0),
   m_backOffDistance(EventSource::s_backoff) {
   std::map<std::string, std::string> pars;
   facilities::Util::keyValueTokenize(params, ", ", pars);
   
   std::string input_file = pars["input_file"];
   ::readLines(input_file, m_inputLines, "#");
   m_interval = 1./std::atof(pars["rate"].c_str());
   if (pars.count("backoff_distance")) {
// NB: Instrument units are in millimeters.
      m_backOffDistance = std::atof(pars["backoff_distance"].c_str());
   }
   m_launchDirection = new FileLaunchDir();
   m_launchPoint = new FileLaunchPoint();
}

auto FileSource::operator() (float xi) -> float {
   (void)(xi);
   return energy(0);
}

auto FileSource::energy(double time) -> double {
   (void)(time);
   if (m_currentLine < m_inputLines.size()) {
      parseCurrentLine();
      m_currentLine++;
   } else {
      // If EOF, use particle properties given in last line of input file.
   }
   return m_energy;
}

void FileSource::parseCurrentLine() {
   std::vector<std::string> tokens;
   facilities::Util::stringTokenize(m_inputLines.at(m_currentLine), 
                                    " \t", tokens);
//   std::cout << m_inputLines.at(m_currentLine) << std::endl;
   m_particleName = tokens.at(0);
   m_energy = std::atof(tokens.at(1).c_str());

// Starting point of incident particle in instrument coordinates.
   double x = std::atof(tokens.at(2).c_str());
   double y = std::atof(tokens.at(3).c_str());
   double z = std::atof(tokens.at(4).c_str());

// Direction cosines of incident particle.
   double cosx = std::atof(tokens.at(5).c_str());
   double cosy = std::atof(tokens.at(6).c_str());
   double cosz = std::atof(tokens.at(7).c_str());

   CLHEP::Hep3Vector dir(cosx, cosy, cosz);
   m_launchDirection->setDir(dir);

   CLHEP::Hep3Vector starting_pt(x, y, z);
   m_launchPoint->setPoint(starting_pt - m_backOffDistance*dir);
}

auto FileSource::particleName() const -> const char * {
   return m_particleName.c_str();
}

auto FileSource::interval(double time) -> double {
   (void)(time);
   return m_interval;
}

auto FileSource::launchDirection() -> LaunchDirection * {
   return m_launchDirection;
}

auto FileSource::launchPoint() -> LaunchPoint * {
   return m_launchPoint;
}

void FileSource::FileLaunchDir::execute(double KE, double time) {
   (void)(KE);
#if 0
   m_glastToGalactic = astro::GPS::instance()->transformGlastToGalactic(time);
#else
   m_glastToGalactic = astro::GPS::instance()->transformToGlast(time, astro::GPS::CELESTIAL).inverse();
#endif
}

auto FileSource::FileLaunchDir::dir() const -> const CLHEP::Hep3Vector & {
   return m_dir;
}

void FileSource::FileLaunchDir::setDir(const CLHEP::Hep3Vector & dir) {
   m_dir = dir;
}

auto FileSource::FileLaunchDir::title() const -> std::string {
   return "FileSource";
}

auto FileSource::FileLaunchDir::skyDirection() const -> astro::SkyDir {
   static CLHEP::Hep3Vector my_dir;
   my_dir = m_glastToGalactic*m_dir;
   return astro::SkyDir(my_dir);
}

auto FileSource::FileLaunchPoint::point() const -> const CLHEP::Hep3Vector & {
   return m_pt;
}

void FileSource::FileLaunchPoint::setPoint(const CLHEP::Hep3Vector & pt) {
   m_pt = pt;
}

auto FileSource::FileLaunchPoint::title() const -> std::string {
   return "FileLaunchPoint";
}
