/** 
* \class ISpectrumFactory
*
* \brief This is an abstract base class for the SpectrumFactory template classes.
* 
* \class RemoteSpectrumFactory
*
* \brief Template class designed to hold method by which to polymorphically instantiate new Spectrum objects
*
* $Header $
*/
#if !defined(AFX_ISPECTRUMFACTORY_H__10BF6F57_E416_4B69_A7CF_220E55281676__INCLUDED_)
#define AFX_ISPECTRUMFACTORY_H__10BF6F57_E416_4B69_A7CF_220E55281676__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000
#include "ISpectrum.h"
//#include "IFluxSvc.h"

#include <typeinfo>
#include <vector>

class ISpectrumFactory
{
public:
    
    /// the only thing this does: make an associated Spectrum object
    virtual ISpectrum* instantiate(const std::string& params)const=0;
    
    /// a dummy to force creation of the spectrum object.
    virtual void addRef()const=0;
    
};


template <class T> class RemoteSpectrumFactory : public ISpectrumFactory 
{
public:
/*    
    RemoteSpectrumFactory(IFluxSvc* svc){
        //Get class name using RTTI:
        std::string classname = typeid(T).name();
        int s = classname.find_first_of("class");
        if( s==0 ) s=6; //found it
        else s =classname.find_first_not_of("0123456789");
        classname = classname.substr(s);
        svc->addFactory(classname, this); 
    }
*/
    /*! return a new Spectrum object
    @param params String from the xml
    @param engine random engine to use
    */

    virtual ISpectrum* instantiate(const std::string& params) const{
        return new T(params);
    }
    
    //! dummy to follow Gaudi model
    virtual void addRef()const{}
};


#endif // !defined(AFX_ISPECTRUMFACTORY_H__10BF6F57_E416_4B69_A7CF_220E55281676__INCLUDED_)
