/** 
* @file SpectrumFactory.h
* @brief declare SpectrumFactory
*
*  $Header$
*/

#ifndef flux_SpectrumFactory_h
#define flux_SpectrumFactory_h

#include "flux/ISpectrumFactory.h"
#include "SpectrumFactoryTable.h"
#include <typeinfo>
#include <vector>

/** 
* \class SpectrumFactory
*
* \brief Template class designed to hold method by which to polymorphically instantiate new Spectrum objects 
* 
*/

template <class T> class SpectrumFactory : public ISpectrumFactory 
{
public:
    
    SpectrumFactory(){
        //Get class name using RTTI:
        m_classname = typeid(T).name();
        int s = m_classname.find_first_of("class");
        if( s==0 ) s=6; //found it
        else s =m_classname.find_first_not_of("0123456789");
        m_classname = m_classname.substr(s);
        SpectrumFactoryTable::instance()->addFactory(m_classname, this); 
    }
    //! return a new Spectrum object
    virtual ISpectrum* instantiate(const std::string& params) const{return new T(params);}
    
    //! dummy to follow Gaudi model
    virtual void addRef()const{}

    std::string name()const{return m_classname;}
private:
    std::string m_classname;
};



#endif
