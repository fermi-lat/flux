// $Header: /nfs/slac/g/glast/ground/cvs/FluxSvc/src/SpectrumFactory.h,v 1.4 2002/07/23 19:00:57 srobinsn Exp $

#if !defined(AFX_SPECTRUMFACTORY_H__211C2F25_9111_44B9_B357_0762789222AF__INCLUDED_)
#define AFX_SPECTRUMFACTORY_H__211C2F25_9111_44B9_B357_0762789222AF__INCLUDED_
/** 
* \class SpectrumFactory
*
* \brief Template class designed to hold method by which to polymorphically instantiate new Spectrum objects 
* 
* $Header $
*/
#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "flux/ISpectrumFactory.h"
#include "SpectrumFactoryTable.h"
#include <typeinfo>
#include <vector>


template <class T> class SpectrumFactory : public ISpectrumFactory 
{
public:
    
    SpectrumFactory(){
        //Get class name using RTTI:
        std::string classname = typeid(T).name();
        int s = classname.find_first_of("class");
        if( s==0 ) s=6; //found it
        else s =classname.find_first_not_of("0123456789");
        classname = classname.substr(s);
        SpectrumFactoryTable::instance()->addFactory(classname, this); 
    }
    //! return a new Spectrum object
    virtual ISpectrum* instantiate(const std::string& params) const{return new T(params);}
    
    //! dummy to follow Gaudi model
    virtual void addRef()const{}
};



#endif // !defined(AFX_SPECTRUMFACTORY_H__211C2F25_9111_44B9_B357_0762789222AF__INCLUDED_)
