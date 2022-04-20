add_library(
  flux STATIC src/CompositeSource.cxx src/Earth.cxx src/EventSource.cxx
  src/FileSource.cxx src/Flux.cxx src/FluxMgr.cxx src/FluxSource.cxx
  src/SimpleSpectrum.cxx src/SourceDirection.cxx src/Spectrum.cxx
  src/SpectrumFactoryTable.cxx src/SurfaceMuons.cxx src/TimeCandle.cxx
  src/VdgGamma.cxx
  src/rootplot/rootAngleHist.cxx src/rootplot/rootAngleHist.h
  src/rootplot/rootEnergyHist.cxx src/rootplot/rootEnergyHist.h
  src/rootplot/rootHist.cxx src/rootplot/rootHist.h src/rootplot/rootplot.cxx
)
target_include_directories(
  flux PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
  $<INSTALL_INTERFACE:>
)
target_link_libraries(
  flux PUBLIC
  astro
  CLHEP::CLHEP
  facilities
  PkgConfig::cfitsio
  xmlBase
)

###### Tests ######
add_executable(test_flux src/test/testMgr.cxx)
target_link_libraries(test_flux PRIVATE flux)

###############################################################
# Installation
###############################################################
install(DIRECTORY flux DESTINATION ${FERMI_INSTALL_INCLUDEDIR})
install(DIRECTORY xml/ DESTINATION ${FERMI_INSTALL_XMLDIR})

install(
  TARGETS flux test_flux
  EXPORT fermiTargets
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
  INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)
