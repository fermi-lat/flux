#ifndef MapSpectrum_H
#define MapSpectrum_H
// $Heading:$
//need fixes:
//1:an assortment of files should be gotten and tested (EGRET skymap)

#include "flux/Spectrum.h"
#include "CLHEP/Random/RandFlat.h"
#include <vector>
#include <string>
#include <map>


//!  Spectrum that reads its differential spectrum from a table
class MapSpectrum : public Spectrum
{
public:
    ///holds information about the intensity and spectral index of a cell
    typedef struct{
        double intensity;
        double index;
    }cellinfo;

    ///Key for indexing the map of cells representing the sky
    class Key{
    public:
        operator int() const{return m_key;}
        Key(double l,double b,double binsize){
            m_key = static_cast<int>( ((360/binsize)*(l+180))+((b+90)/binsize) );}
        void init(double l,double b,double binsize){
            m_key = static_cast<int>( ((360/binsize)*(l+180))+((b+90)/binsize) );}
    private:
        int m_key;
    };

    class RowKey{
    public:
        operator int() const{return m_key;}
        RowKey(double b,double binsize){
            m_key = static_cast<int>( ((180/binsize)*b) );}
    private:
        int m_key;
    };

    /// params is the filename to read
    MapSpectrum(const std::string& params);
    /// default ctor for factoryTable
    MapSpectrum();


    /// return total flux 
    virtual double flux(double time) const;

    ///defaulted interval() to -1.
    virtual double interval(double time);

    /// sample a single particle energy from the spectrum
    virtual float operator() (float);

    ///this returns a galactic direction, in the form (l,b)
    std::pair<double,double> dir(double energy);

    ///this returns a galactic direction, in the form (l,b)
    std::pair<double,double> findDir(double energy);

    virtual std::string title() const;
    virtual const char * particleName() const;
    inline  const char * nameOf() const {return "MapSpectrum";}
    //   use default destructor, copy constructor, and assignment op.

    /// returns the solid angular size of a bin of angular dinension
    /// "m_binSize", located at colatitude b
    double sizeOfBin(double b);

    ///sets the net flux
    void setNetFlux();   
    float parseParamList(std::string input, int index);

    ///stuff required by IService, should be vestigial
    virtual std::pair<float,float> dir(float energy)const;
    virtual double energy(double time=0);


private:
    ///the map of cells covering the sky
    std::map<Key,cellinfo> m_catalog;
    std::vector<cellinfo> m_testCatalog;

    ///the map of scaled intensities for each row of the map.
    std::map<RowKey,double> m_rowCatalog;


    double m_flux;   // current flux (set when direction changes)
    double m_binSize; // bin size of the sky tessellation (this should come from the FITS file)

    double m_currentl,m_currentb; //current location on sky of incoming photon

    double m_index;  //current power-law index for the energy (set when direction changes)
    double m_netFlux; //sum over the entire sky's flux
    //float m_E0;  //energy base
    float m_mapE0; //lowest energy (contained in map).
    float m_extrapE0; // lowest energy (to be used - extrapolated).
    float m_EMax; // highest energy (cutoff for sending photons)
    float m_scaleFactor; //factor by which to scale all the flux bins(hence,the whole flux)


    //default spectral index
    double m_defaultIndex;

    std::string initialization_document;
    std::string m_particle_name;
};



#endif // MapSpectrum_H
