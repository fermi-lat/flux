def generate(env, **kw):
    env.Tool('addLibrary', library = ['flux'], package = 'flux')
    env.Tool('xmlBaseLib')
    env.Tool('astroLib')
    env.Tool('addLibrary', library = env['clhepLibs'])
    env.Tool('addLibrary', library = env['cfitsioLibs'])

def exists(env):
    return 1
