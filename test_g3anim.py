import numpy as np
from g3suite.g3anim import *
from g3suite import g3read_units as g3u,  g3maps, g3matcha, g3anim
import g3suite.g3read as g3read
g3anim.debug = True
g3matcha.cache_filename = 'cache.pk2'
g3matcha.cache_type = 'shelve' #pickle'
g3matcha.cache_from_filename_only = True
g3matcha.size_limit = 5200
g3matcha.debug=True
g3maps.debug = True
g3u.debug=True
g3.debug = True

snap_base = './'

snap_indexes = ['%03d'%i for i in range(20,93)]


snaps = ['snapdir_%s/snap_%s'%(s,s) for s in snap_indexes]
groups = ['groups_%s/sub_%s'%(s,s) for s in snap_indexes]
fofs = ['groups_%s/group_tab_%s'%(s,s) for s in snap_indexes]





units = g3u.get_units(g3.get_one_file(snap_base+snaps[0]),  time=1.)
ureg = units.get_u()

anim = Animation(fps=200)

camera = anim.add_character(CameraGadgetParticlesFollow(units = units, ureg=ureg, file_format='movie/%04d.png', dmo=True))
keyframes = anim.add_character(CharacterGadgetParticles(units = units, ureg=ureg, dmo=True))
size = 2e3*ureg.kpc
t0=None

@g3matcha.memoize
def read_fofs(filename, boxsize=None, use_cache = False):
    return g3read.read_fofs(filename, boxsize = boxsize)
boxsize = None
for i in range(0,len(snaps)):
    if boxsize==None:
        f = g3read.GadgetFile(snap_base+snaps[i]+'.0')
        poses = f.read_new('POS ',-1)
        boxsize = f.header.BoxSize
        print('pos min',np.min(poses), 'pos max', np.max(poses), 'boxsize',boxsize)
    print('base:',  snap_base+groups[i])
    d = read_fofs(snap_base + fofs[i], boxsize = boxsize, use_cache = True)
    #print(d)
    ds = g3matcha.numpy_to_dict(d, blocks=['GPOS','GVEL','len'])
    #print(ds)
    
    ds.sort(
        key = lambda x :-x['len'])
    #print(ds)
    first_fof = ds[0]
    print(first_fof)
    halo_pos = first_fof['GPOS']
    halo_vel = first_fof['GVEL']
    print(halo_pos, halo_vel)

    kf = keyframes.add_snapshot(snap_base+snaps[i])
    t = kf['t']
    if t0 is None:
        t0 = t
    print('t',t)
    camera.add_keyframe(t,{"GPOS":halo_pos*ureg.glength, "GVEL":halo_vel*ureg.cvelocity, "SIZE":size})

print('red!')



anim.anim(camera, t, t_start=t0, i_start=537)
