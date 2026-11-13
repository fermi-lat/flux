//
//////////////////////////////////////////////////////////////////////

#include "flux/SpectrumFactoryTable.h"

#include <list>
#include <string>

SpectrumFactoryTable* SpectrumFactoryTable::s_instance = nullptr;

auto SpectrumFactoryTable::instantiate(const std::string& name, 
                                             const std::string& params) const -> ISpectrum* 
{
    auto returnloc = find(name);
    return returnloc==end() ? nullptr: (*returnloc).second->instantiate(params);
}
auto SpectrumFactoryTable::instantiate(const std::string& name) const -> ISpectrum* 
{
    auto returnloc = find(name);
    std::string params; // dummy 
    return returnloc==end() ? nullptr: (*returnloc).second->instantiate(params);
}


auto SpectrumFactoryTable::spectrumList()const -> std::list<std::string>
{
    
    std::list<std::string>  outstr;
    for( auto tableiter = begin(); tableiter!=end() ; ++tableiter){
        outstr.push_back(tableiter->first);
    }
    return outstr;
    
}
