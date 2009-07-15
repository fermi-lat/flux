# -*- python -*-
# $Id: SConscript,v 1.8 2009/07/14 23:23:44 jrb Exp $
# Authors: Toby Burnett <tburnett@u.washington.edu>
# Version: flux-08-40-03
Import('baseEnv')
Import('listFiles')
progEnv = baseEnv.Clone()
libEnv = baseEnv.Clone()

libEnv.Tool('fluxLib', depsOnly = 1)
fluxLib = libEnv.StaticLibrary('flux', listFiles(['src/*.cxx','src/rootplot/*.cxx']))

progEnv.Tool('fluxLib')
test_fluxBin = progEnv.Program('test_flux',[ 'src/test/testMgr.cxx'])
rootplotsBin = progEnv.Program('rootplots',[ 'src/rootplot/test/roottest.cxx'])

#progEnv.Tool('registerObjects', package = 'flux', libraries = [fluxLib], binaries = [rootplotsBin], testApps = [test_fluxBin], includes = listFiles(['flux/*.h']),
#             xml = listFiles(['xml/*'], recursive = True))

progEnv.Tool('registerTargets', package = 'flux',
             staticLibraryCxts = [[fluxLib,libEnv]],
             binaryCxts = [[rootplotsBin, progEnv]],
             testAppCxts = [[test_fluxBin, progEnv]],
             includes = listFiles(['flux/*.h']),
             xml = listFiles(['xml/*'], recursive = True))




