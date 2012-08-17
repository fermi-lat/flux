#$Id: fluxLib.py,v 1.2 2008/02/26 03:23:33 glastrm Exp $
def generate(env, **kw):
    if not kw.get('depsOnly',0):
        env.Tool('addLibrary', library = ['flux'])
	if env['PLATFORM'] == 'win32' and env.get('CONTAINERNAME','') == 'GlastRelease':
	    env.Tool('findPkgPath', package = 'flux') 
	    env.Tool('findPkgPath', package = 'astro') 
	    env.Tool('findPkgPath', package = 'facilities') 
    env.Tool('xmlBaseLib')
    env.Tool('astroLib')
    env.Tool('addLibrary', library = env['clhepLibs'])
    env.Tool('addLibrary', library = env['cfitsioLibs'])
    if kw.get('incsOnly', 0) == 1: 
        env.Tool('findPkgPath', package = 'facilities') 
        env.Tool('findPkgPath', package = 'astro')
        env.Tool('findPkgPath', package = 'xmlBase')

def exists(env):
    return 1
