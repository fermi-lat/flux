add_library(flux STATIC src/CompositeSource.cxx src/Earth.cxx src/EventSource.cxx
      src/FileSource.cxx src/Flux.cxx src/FluxMgr.cxx src/FluxSource.cxx
      src/SimpleSpectrum.cxx src/SourceDirection.cxx src/Spectrum.cxx
      src/SpectrumFactoryTable.cxx src/SurfaceMuons.cxx src/TimeCandle.cxx
      src/VdgGamma.cxx
      src/rootplot/rootAngleHist.cxx src/rootplot/rootAngleHist.h
      src/rootplot/rootEnergyHist.cxx src/rootplot/rootEnergyHist.h
      src/rootplot/rootHist.cxx src/rootplot/rootHist.h src/rootplot/rootplot.cxx
      )
target_include_directories(flux PUBLIC
                           $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
                           $<INSTALL_INTERFACE:>)
target_link_libraries(flux PUBLIC
    astro
    CLHEP::CLHEP
    facilities
    PkgConfig::cfitsio
    xmlBase)

# target_compile_definitions(flux PUBLIC BUILD_WITHOUT_ROOT)

# add_library(EarthPhenom STATIC
#     EarthPhenom/src/EarthPhenomInner.cxx
#     EarthPhenom/src/EarthPhenomLimb.cxx)
# target_include_directories(EarthPhenom PUBLIC
#                            $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/EarthPhenom>
#                            $<INSTALL_INTERFACE:>)
# target_link_libraries(EarthPhenom PUBLIC flux astro facilities)

# add_subdirectory(eblAtten)
# add_subdirectory(genericSources)
# add_subdirectory(GRB)
# add_subdirectory(GRBobs)
# add_subdirectory(GRBtemplate)
# add_subdirectory(microQuasar)
# add_subdirectory(Pulsar)
# add_subdirectory(pySources)
# add_subdirectory(SpectObj)
