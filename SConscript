# -*- python -*-
# $Id: SConscript,v 1.15 2010/06/12 22:43:36 jrb Exp $
# Authors: Toby Burnett <tburnett@u.washington.edu>
# Version: flux-08-40-08
Import('baseEnv')
Import('listFiles')
progEnv = baseEnv.Clone()
libEnv = baseEnv.Clone()

fluxLib = libEnv.StaticLibrary('flux',
                               listFiles(['src/*.cxx','src/rootplot/*.cxx']))

progEnv.Tool('fluxLib')
test_fluxBin = progEnv.Program('test_flux',[ 'src/test/testMgr.cxx'])
rootplotsBin = progEnv.Program('rootplots',[ 'src/rootplot/test/roottest.cxx'])


progEnv.Tool('registerTargets', package = 'flux',
             staticLibraryCxts = [[fluxLib,libEnv]],
             binaryCxts = [[rootplotsBin, progEnv]],
             testAppCxts = [[test_fluxBin, progEnv]],
             includes = listFiles(['flux/*.h']),
             xml = listFiles(['xml/*'], recursive = True))




