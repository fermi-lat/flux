# -*- python -*-
# $Id: SConscript,v 1.3 2008/03/21 16:52:34 glastrm Exp $
# Authors: Toby Burnett <tburnett@u.washington.edu>
# Version: flux-08-40-02
Import('baseEnv')
Import('listFiles')
progEnv = baseEnv.Clone()
libEnv = baseEnv.Clone()

libEnv.Tool('fluxLib', depsOnly = 1)
fluxLib = libEnv.StaticLibrary('flux', listFiles(['src/*.cxx','src/rootplot/*.cxx']))

progEnv.Tool('fluxLib')
test_fluxBin = progEnv.Program('test_flux', 'src/test/testMgr.cxx')
rootplotsBin = progEnv.Program('rootplots', 'src/rootplot/test/roottest.cxx')

progEnv.Tool('registerObjects', package = 'flux', libraries = [fluxLib], binaries = [rootplotsBin], testApps = [test_fluxBin], includes = listFiles(['flux/*.h']),
             xml = listFiles(['xml/*'], recursive = True))
