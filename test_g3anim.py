import g3anim
g3anim.debug = True
from g3anim import *
import g3read_units as g3u, g3read as g3, numpy as np
import g3maps
import g3matcha
g3matcha.cache_filename = 'cache.pk2'
g3matcha.cache_from_filename_only = True
g3matcha.size_limit = 5200
g3matcha.debug=True
g3maps.debug = True
g3u.debug=True
g3.debug = True

snap_base = '/HydroSims/Magneticum/Box4/uhr_test/'
snaps = ['snapdir_012/snap_012', 'snapdir_016/snap_016', 'snapdir_020/snap_020', 'snapdir_024/snap_024', 'snapdir_028/snap_028', 'snapdir_032/snap_032', 'snapdir_036/snap_036', 'snapdir_040/snap_040', 'snapdir_044/snap_044', 'snapdir_048/snap_048', 'snapdir_052/snap_052', 'snapdir_058/snap_058', 'snapdir_060/snap_060', 'snapdir_064/snap_064', 'snapdir_068/snap_068', 'snapdir_072/snap_072', 'snapdir_076/snap_076', 'snapdir_080/snap_080', 'snapdir_084/snap_084', 'snapdir_088/snap_088', 'snapdir_092/snap_092', 'snapdir_096/snap_096', 'snapdir_100/snap_100', 'snapdir_102/snap_102', 'snapdir_104/snap_104', 'snapdir_106/snap_106', 'snapdir_108/snap_108', 'snapdir_112/snap_112', 'snapdir_116/snap_116', 'snapdir_120/snap_120', 'snapdir_124/snap_124', 'snapdir_128/snap_128', 'snapdir_132/snap_132', 'snapdir_136/snap_136']

groups = [snap.replace('snapdir','groups').replace('snap','sub') for snap in snaps]


units = g3u.get_units(g3.get_one_file(snap_base+snaps[0]),  time=1.)
ureg = units.get_u()

anim = Animation(fps=50)

camera = anim.add_character(CameraGadgetParticlesFollow(units = units, ureg=ureg, file_format='frames_2/%04d.png'))
keyframes = anim.add_character(CharacterGadgetParticles(units = units, ureg=ureg))
size = .5e3*ureg.kpc
t0=None
print('base:',  snap_base+groups[0])
g = g3matcha.yield_haloes(snap_base+groups[0], ihalo_start=0, ihalo_end=2, use_cache=True)
next(g)
first_halo = next(g)
ihalo = first_halo['ihalo']
halo_pos = first_halo_pos = first_halo['GPOS']
kf = keyframes.add_snapshot(snap_base+snaps[0])
t = kf['t']
if t0 is None:
    t0 = t
#data = keyframes.memo_read(kf['path'], first_halo_pos, size, self.blocks, -1, use_cache = True)

subhalo_data_vel = next(g3matcha.yield_subhaloes(snap_base + groups[0],  ihalo,   use_cache = True, blocks=('SVEL',), halo_goff= first_halo['GOFF']))['SVEL']
halo_vel = subhalo_data_vel
halo1_pos = None
halo1_vel = None
t1 = None
print('#halo pos ',halo_pos, ' halo vel ',halo_vel)
camera.add_keyframe(t,{"GPOS":halo_pos*ureg.kpc, "GVEL":halo_vel*ureg.cvelocity, "SIZE":size})

for i in range(1,len(snaps)):
#for i in range(1,3):
    halo = next(g3matcha.yield_matches(snap_base+snaps[i-1], halo_pos, size,
                  snap_base + groups[i], snap_base + groups[i],
                  None, None, 1e2, None, size*10., None, use_cache = True))
    halo_pos = halo['GPOS']
    ihalo = halo['ihalo']
    subhalo_data_vel = next(g3matcha.yield_subhaloes(snap_base + groups[i], ihalo, use_cache = True, blocks=('SVEL',), halo_goff = halo['GOFF']))['SVEL']
    halo_vel = subhalo_data_vel
    kf = keyframes.add_snapshot(snap_base+snaps[i])
    t = kf['t']
    
    halo1_pos = halo_pos
    halo1_vel = halo_vel
    #data = keyframes.memo_read(kf['path'], first_halo_pos, size, self.blocks, -1, use_cache = True)
    a = keyframes.get_scale_factor(t)
    if a<0.2 or a>0.4:
        camera.add_keyframe(t,{"GPOS":halo_pos*ureg.kpc, "GVEL":halo_vel*ureg.cvelocity, "SIZE":size})
    
#camera.add_keyframe(t0,{"GPOS":halo_pos*ureg.kpc, "GVEL":np.array([0.,0.,0.])*ureg.cvelocity, "SIZE":size})
#





#1/0
anim.anim(camera, t, t_start=4., i_start=0)
