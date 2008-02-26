#$Id$
def generate(env, **kw):
    if not kw.get('depsOnly',0):
        env.Tool('addLibrary', library = ['flux'])
    env.Tool('xmlBaseLib')
    env.Tool('astroLib')
    env.Tool('addLibrary', library = env['clhepLibs'])
    env.Tool('addLibrary', library = env['cfitsioLibs'])

def exists(env):
    return 1
