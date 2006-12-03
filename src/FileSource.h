/**
 * @file FileSource.h
 * @brief This source reads the incident particle properties -- 
 * particle type, 4-momentum, trajectory -- from a file.
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/flux/src/FileSource.h,v 1.4 2006/03/21 01:28:56 usher Exp $
 */

#ifndef _flux_FileSource_h
#define _flux_FileSource_h

#include <string>
#include <utility>
#include <vector>

#include "flux/LaunchDirection.h"
#include "flux/LaunchPoint.h"
#include "flux/Spectrum.h"

/**
 * @class FileSource
 * @brief This source reads the incident particle properties -- 
 * particle type, 4-momentum, trajectory -- from a file.
 *
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/flux/src/FileSource.h,v 1.4 2006/03/21 01:28:56 usher Exp $
 */

class FileSource : public Spectrum {

public:

   FileSource(const std::string & params);

   virtual ~FileSource() {
      delete m_launchDirection;
      delete m_launchPoint;
   }

   /// @return Particle energy (MeV)
   /// @param xi Uniform random deviate on the unit interval.
   virtual float operator() (float xi);

   /// @return Particle type
   virtual const char * particleName() const;

   /// @return Title describing the spectrum.
   virtual std::string title() const {
      return "FileSource";
   }

   /// @return Interval to the next event (seconds)
   virtual double interval(double time);

   /// @return Photon energy (MeV).
   virtual double energy(double time);

   /// @return Pointer to the LaunchDirection data member.
   virtual LaunchDirection * launchDirection();

   /// @return Pointer to the LaunchPoint data member.
   virtual LaunchPoint * launchPoint();

protected:

   /// Disable these virtual functions since they are not used by
   /// this source.
   virtual double flux(double) const {return 0;}
   virtual double solidAngle() const {return 0;}
   virtual std::pair<double, double> dir(double) {
      return std::make_pair(0, 0);
   }

private:
   
#ifndef SWIG
   class FileLaunchDir : public LaunchDirection {
   public:
      FileLaunchDir() {}
      virtual ~FileLaunchDir() {}
      virtual void execute(double KE, double time);

      /// @return The particle direction in instrument coordinates.
      virtual const CLHEP::Hep3Vector & dir() const;

      virtual std::string title() const;
      virtual astro::SkyDir  skyDirection() const;
      virtual double zenithCosine() const {
         return 1.;
      }

      void setDir(const CLHEP::Hep3Vector & dir);
   private:
      CLHEP::Hep3Vector m_dir;
      CLHEP::HepRotation m_glastToGalactic;
   } * m_launchDirection;

   class FileLaunchPoint : public LaunchPoint {
   public:
      FileLaunchPoint() {}
      virtual ~FileLaunchPoint() {}
      virtual const CLHEP::Hep3Vector & point() const;
      virtual std::string title() const;

      void setPoint(const CLHEP::Hep3Vector & pt);
   private:
      CLHEP::Hep3Vector m_pt;
   } * m_launchPoint;
#endif // SWIG

   std::vector<std::string> m_inputLines;
   unsigned int m_currentLine;
   double m_interval;

   double m_backOffDistance;

   double m_energy;
   std::string m_particleName;

   void parseCurrentLine();

};

#endif // _flux_FileSource_h
