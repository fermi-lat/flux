//
//////////////////////////////////////////////////////////////////////

#include "flux/SpectrumFactoryTable.h"

#include <list>
#include <string>

SpectrumFactoryTable* SpectrumFactoryTable::s_instance = 0;

auto SpectrumFactoryTable::instantiate(const std::string& name, 
                                             const std::string& params) const -> ISpectrum* 
{
    const_iterator returnloc = find(name);
    return returnloc==end() ? 0: (*returnloc).second->instantiate(params);
}
auto SpectrumFactoryTable::instantiate(const std::string& name) const -> ISpectrum* 
{
    const_iterator returnloc = find(name);
    std::string params; // dummy 
    return returnloc==end() ? 0: (*returnloc).second->instantiate(params);
}


auto SpectrumFactoryTable::spectrumList()const -> std::list<std::string>
{
    
    std::list<std::string>  outstr;
    for( const_iterator tableiter = begin(); tableiter!=end() ; ++tableiter){
        outstr.push_back(tableiter->first);
    }
    return outstr;
    
}
