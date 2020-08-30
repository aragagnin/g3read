from g3anim import *
import g3read_units as g3u, g3read as g3, numpy as np
import g3maps
import g3matcha

g3matcha.debug=True
g3maps.debug = True
g3u.debug=True
g3.debug = True

snap_base = '/HydroSims/Magneticum/Box4/uhr_test/'
snaps = ['snapdir_012/snap_012', 'snapdir_016/snap_016', 'snapdir_020/snap_020', 'snapdir_024/snap_024', 'snapdir_028/snap_028', 'snapdir_032/snap_032', 'snapdir_036/snap_036', 'snapdir_040/snap_040', 'snapdir_044/snap_044', 'snapdir_048/snap_048', 'snapdir_052/snap_052', 'snapdir_058/snap_058', 'snapdir_060/snap_060', 'snapdir_064/snap_064', 'snapdir_068/snap_068', 'snapdir_072/snap_072', 'snapdir_076/snap_076', 'snapdir_080/snap_080', 'snapdir_084/snap_084', 'snapdir_088/snap_088', 'snapdir_092/snap_092', 'snapdir_096/snap_096', 'snapdir_100/snap_100', 'snapdir_102/snap_102', 'snapdir_104/snap_104', 'snapdir_106/snap_106', 'snapdir_108/snap_108', 'snapdir_112/snap_112', 'snapdir_116/snap_116', 'snapdir_120/snap_120', 'snapdir_124/snap_124', 'snapdir_128/snap_128', 'snapdir_132/snap_132', 'snapdir_136/snap_136']

units = g3u.get_units(g3.get_one_file(snap_base+snaps[0]))
ureg = units.get_u()

center_cu = np.array([32431.984375, 33903.015625, 26402.964844])*ureg.glength
radius_cu = 1.e5*ureg.glength

anim = Animation()

camera = anim.add_character(CameraGadgetParticlesStill(center_cu, radius_cu))
keyframes = anim.add_character(CharacterGadgetParticles(units = units, ureg=ureg))
for snap in snaps:
    keyframes.add_snapshot(snap_base+snap)

anim.anim(10.)
