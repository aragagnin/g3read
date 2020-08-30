import bisect, math, g3read_units, g3maps, g3read_units as g3u, matplotlib.pyplot as plt, numba, g3read as g3, g3matcha

class KeyFrames(list):
    def insert(self, t, values):
        new_entry = dict(values)
        new_entry["t"] = t
        keys = [d['t'] for d in self]
        return super(KeyFrames, self).insert(bisect.bisect_left(keys, new_entry['t']), new_entry)

class Character(object):
    def __init__(self, props = None):
        self._props = props or []
        self._keyframes = KeyFrames()
        self._interp =  Character.default_interp
        self._interp_values = Character.default_interp_values
        self._is_camera = False
        self._dict_to_frame = self.default_dict_to_frame
    def default_dict_to_frame(self, d):
        return d
    
    def default_interp_values(v1,v2,t1,t2,t):
        delta_t = t2-t1
        factor = (t-t1)/delta_t
        v = v1 + ( v2 -  v1)*factor
        #print('v1',v1,'v1',v2,'t1',t1,'t2',t2,'t',t,'dt',delta_t,'fac',factor,'v',v)
        return v
    
    def default_interp(kf1,kf2,t, _interp_values):
        if kf1 is not None:
            t1 = kf1['t']
        if kf2 is not None:
            t2 = kf2['t']
        if kf1 is not None:
            keys =  kf1.keys()
        elif kf2 is not None:
            keys =  kf2.keys()
        else:
            return None
        f = {}
        for key in keys:
            #print(key, kf1, kf2)
            if kf2 is  None or key not in kf2 or kf2[key]==None:
                f[key] = kf1[key]
            elif kf1 is  None or key not in kf1 or kf1[key]==None:
                f[key] = kf2[key]
            else:
                f[key] = _interp_values(kf1[key], kf2[key], t1, t2, t)
        return f
    
    def add_keyframe(self, t, values):
        self._keyframes.insert(t, values)

    def render_start(self):
        pass
    def get_keyframes_nearby(self, i, fps):
        kf_prev = None
        kf_next = None
        kf_this = None
        for kf in self._keyframes:
            _i = math.floor(kf['t']*fps)
            rkf = self._dict_to_frame(kf)
            if _i<i:
                kf_prev = rkf
            elif _i==i:
                kf_this = rkf
            elif kf_next is None:
                kf_next = rkf
                break
        return kf_prev, kf_this, kf_next
            

        
class Camera(Character):
    def __init__(self, props):
        super(Camera, self).__init__(props)
        self._is_camera = True
        self._render = Camera.default_render
        
    def default_render(characters, i, t):
        res = '['
        for character in characters:
            res+= str(character)+','
        print ('t=%f i=%d '%(t,i)+res+']')
    def default_interp_values(v1,v2,t1,t2,t):        
class Animation(object):
    def __init__(self, fps):
        self._characters = []
        self._fps = fps
    def add_character(self, c=None, props=None):
        if c is None:
            c = Character(props)
        self._characters.append(c)
        return c
    
    def anim(self, camera, t=None):
        i = -1
        fps = self._fps
        camera.render_start()
        while True:
            i+=1
            my_t = i/fps
            if t is not None and my_t>t:
                break
            delta_t = 1/fps
            frame = []
            for character in self._characters:
                kf1, kf, kf2 =  character.get_keyframes_nearby(i, fps)
                if kf is not None:
                    frame.append(kf)
                else:
                    f = character._interp(kf1, kf2, my_t, character._interp_values)
                    if f is None:
                        return
                    frame.append(f)
            
            camera._render(frame, i, my_t)
        camera.render_end()

@numba.jit(nopython=True)
def interp_0_1(x0p,  x1, x1p,  t):
    D = 0
    C = x0p
    B = 3.*x1 - x1p - 2.*C - 3*D
    A = x1 - B - C - D
    return A*t*t*t + B*t*t  + C*t + D

@numba.jit(nopython=True)
def interp(t0, x0, x0p, t1, x1, x1p, t):
    delta_x = x1 - x0
    delta_t = t1 - t0
    return x0 + interp_0_1(x0p*delta_t, x1-x0, x1p*delta_t,  (t - t0)/delta_t)


class CameraGadgetParticles(Camera):
    def __init__(self):
        super(Camera, self).__init__(props = ["CPOS","OPOS"])

    def default_render(frames, i, t):
        f,a  = plt.subplot()
        for frame in frames:
            a.scatter(frame['POS '][:,0], frame['POS '][:,0])
        f.save_fig('frame_%4d'%(i))


class CharacterGadgetParticles(Character):
    def __init__(self, blocks=None, time_factor=10.):
        blocks_default = ["POS ","VEL ","MASS","ID  ","RHO "]
        super(Character, self).__init__()
        self.blocks = list(set(blocks + blocks_default)) if blocks is not None else blocks_default #we get uniq values
        self.prop = "MASS"
        self.time_factor = time_factor
        self._interp_values = CharacterGadgetParticles.default_interp_values
        self.om0 = 1.
        self.q0 = 1.
        self.h0 = .7 #TODO: read hubble constant h0 from simulation data
        self.ureg=  None
    def get_time(self, f):
        return f.header.time * self.time_factor
    
    def add_snapshot(self, snap_path):
        
        f = g3.GadgetFile(get_one_file(snap_path))
        if self.ureg is None : #we still ahvent set units
            self.units = g3u.Units();
            self.units.set_context(f.header.HubbleParam, 1., debug=True) #we will fo stuff in comoving so a=1.
        self._keyframes.insert(self.get_time(f), {"path":snap_path})

    # we cache reading of data
    @g3matcha.memoize
    def memo_read(path, pos, size, blocks, ptypes):
        return g3u.read_particles_in_box(path, pos, size, blocks, ptypes, units=self.units)
    
    def default_dict_to_frame(self, d, camera):
        data = self.memo_read(d['path'], camera.gpos_cu, camera.max_xyz_cu, self.blocks, -1)
        ids1 = data['ID  ']
        ids1_argsort = np.argsort(ids1)
        for k in data:
            data[k] = data[k][ids1_argsort]
        return data
    
    @numba.jit(nopython=True)    
    def default_interp_values(t1, ids1, pos1, vel1, t2, ids2, pos2, vel2, t3, ids3, pos3):
        N3 = len(ids3)
        ids1_count = 0
        ids2_count = 0
        for i in range(N3):
            id3 = ids3[i]
            while ids1[ids1_count]<id3:
                ids1_count +=1
            while ids2[ids2_count]<id3:
                ids2_count +=1
            id1 = ids1[ids1_count]
            id2 = ids2[ids2_count]
            if id1 == id3 and id2==id3:
                for j in range(3):
                    pos3[i][j] = interp(t1, pos1[id1i_count,j], vel1[ids1_count,j], t2, pos2[ids2_count,j], vel2[ids2_count,j], t3)
            elif id1 == id3:
                for j in range(3):
                    pos3[i][j] = pos1[ids1_count,j] + (t3-t1)*vel1[ids1_count,j]
            elif id2 == id3:
                for j in range(3):
                    pos3[i][j] = pos2[ids2_count,j] + (t3-t2)*vel2[ids2_count,j]

    def scale_factor_to_time(self, a):
        "formula to convert scale factor in Gyr stolen from the cosmological calculator https://home.fnal.gov/~gnedin/cc/"
        x = (1.00001-self.om0)/self.om0*a*a*a
        w = 651.3333/self.h0/np.sqrt(1.00001-self.om0)*self.q0*np.log(np.sqrt(x)+np.sqrt(x+1));
        return w*self.units.yr*1e9
    def default_interp(kf1,kf2,t3, _interp_values):
        t1 = self.scale_factor_to_time(kf1['t']).to('s').magnitude
        t2 = self.scale_factor_to_time(kf2['t']).to('s').magnitude
        ids1 = kf1['ID  ']
        ids2 = kf2['ID  ']
        pos1 = kf1['POS '].to('kpc').magnitude
        pos2 = kf2['POS '].to('kpc').magnitude
        vel1 = kf1['VEL '].to('kpc/s').magnitude
        vel2 = kf2['VEL  '].to('kpc/s').magnitude
        ids3 = np.sort(np.uniq(np.sort(np.concate(ids1, ids2))))
        N3 = len(ids3)
        pos3 = np.zeros(N3*3).reshape(N3,3)
        default_interp_values(ids1, pos1, vel1, ids2, pos2, vel2, ids3, pos3, t3)
        pos3*=ureg.kpc
        return {"POS ": pos3}
    def default_render(characters, i, t):
        print('we!')


class SnapshotsAnimation(Animation):
    def __init__(self, length = 10., fps=50,  snap_to_seconds = scalefactor_to_seconds):
        super(Animation, self).__init__(fps=fps)
        self.length = length
        self.snap_to_seconds = snap_to_seconds
        
    def add_character(self, c=None, props=None):
        if c is None:
            c = Character(props)
        self._characters.append(c)
        return c

    def add_character(self, props):
        c = Character(props)
        self._characters.append(c)
        return c

    def anim(self, t=None):
        i = -1
        fps = self._fps
        while True:
            i+=1
            my_t = i/fps
            if t is not None and my_t>t:
                break
            delta_t = 1/fps
            frame = []
            for character in self._characters:
                kf =  character.get_keyframe_or_none(i, fps)
                if kf is not None:
                    #print('KF!', kf)
                    frame.append(kf)
                else:
                    kf1 = character.get_prev_keyframe(i, fps)
                    kf2 = character.get_next_keyframe(i, fps)
                    #print('kf1,kf2', kf1, kf2)
                    f = character._interp(kf1, kf2, my_t, character._interp_values)
                    if f is None:
                        return
                    frame.append(f)
            
            yield self._render(frame, i, my_t)
