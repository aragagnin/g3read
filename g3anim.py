import bisect, math, matplotlib.pyplot as plt, numba,  numpy as np
from . import g3maps, g3read_units as g3u, g3read as g3, g3matcha

debug = False

class KeyFrames(list):
    def insert(self, t, values):
        new_entry = dict(values)
        new_entry["t"] = t
        keys = [d['t'] for d in self]
        super(KeyFrames, self).insert(bisect.bisect_left(keys, new_entry['t']), new_entry)
        return new_entry
class Character(object):
    def __init__(self, props = None):
        self._props = props or []
        self.keyframes = KeyFrames()
        self.cache  = None
        self._is_camera = False
        self._dict_to_frame = self.default_dict_to_frame
    def default_dict_to_frame(self, d, camera):
        d['is_keyframe']=False
        return d
    def interp_values(self, v1,v2,t1,t2,t):
        delta_t = t2-t1
        factor = (t-t1)/delta_t
        v = v1 + ( v2 -  v1)*factor
        #print('v1',v1,'v1',v2,'t1',t1,'t2',t2,'t',t,'dt',delta_t,'fac',factor,'v',v)
        return v
    
    def interp(self, kf1,kf2,t):
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
                f[key] = self.interp_values(kf1[key], kf2[key], t1, t2, t)
        return f
    
    def add_keyframe(self, t, values):
        return self.keyframes.insert(t, values)

    def render_start(self):
        pass
    def render_end(self):
        pass
    def get_frame(self, i, camera):
        if self.cache is not None and i in self.cache:
            return self.cache[i]
        
        kf1, kf, kf2 =  self.get_keyframes_nearby(i,  camera)
        print('[%c][%c][%c]'%(
            'x' if kf1 is not None else ' ',
            'x' if kf is not None else ' ',
            'x' if kf2 is not None else ' '))
        if kf is not None:
            f  = kf
        else:
            t = kf1['t'] + (kf2['t'] - kf1['t'])*(i-kf1['i'])/(kf2['i']-kf1['i'])
            print('# call self.interp()')
            f = self.interp(kf1, kf2, t)
            print('# called self.interp()')
        if self.cache is not None:
            self.cache[i] = f
        return f
    def get_keyframes_nearby(self, i,  camera):
        kf_prev = None
        kf_next = None
        kf_this = None

        i_prev = None
        i_this = None
        i_next = None
        for kf in self.keyframes:
            _i = kf['i']
            if _i<i:
                i_prev = _i
                kf_prev = kf
            elif _i == i:
                i_this = _i
                kf_this = kf
            elif i_next is None:
                i_next = _i
                kf_next = kf
                break
        print('# dicts to kv_prev, kv_this, kv_next - has camera?', camera is not None)
        if kf_prev is not None:        
            camera_d = None

            if camera is not None:
                print('# prev get frame')
                camera_d =  camera.get_frame(i_prev, None)
            print('# prev dict to  frame')
            kf_prev = self._dict_to_frame(kf_prev, camera_d)
            kf_prev['is_keyframe']=True
        if kf_this is not None:        
            camera_d = None
            if camera is not None:
                print('# this get frame')
                camera_d =  camera.get_frame(i_this, None)
            print('# this dict to  frame')
            kf_this = self._dict_to_frame(kf_this, camera_d)
            kf_this['is_keyframe']=True
        if kf_next is not None:        
            camera_d = None
            if camera is not None:
                print('# bext get frame')
                camera_d =  camera.get_frame(i_next, None)
            print('# next dict to  frame')
            kf_next = self._dict_to_frame(kf_next, camera_d)
            kf_next['is_keyframe']=True
        return kf_prev, kf_this, kf_next

            

        
class Camera(Character):
    def __init__(self, props = None):
        super(Camera, self).__init__(props=props)
        self._is_camera = True

        
    def render(self,characters, i, t):
        res = '['
        for character in characters:
            res+= str(character)+','
        print ('t=%f i=%d '%(t,i)+res+']')

        
class Animation(object):
    def __init__(self, fps=50):
        self._characters = []
        self._fps = fps
    def add_character(self, c=None, props=None):
        if c is None:
            c = Character(props=props)
        self._characters.append(c)
        return c
    
    def anim(self, camera, t=None, t_start = 0., i_start=0):
        i = -1
        fps = self._fps
        
        for character in self._characters:
            for kf in character.keyframes:
                kf['i'] = math.floor(kf['t']*self._fps)
                
        camera.render_start()
        while True:
            i+=1
            my_t = i/fps
            if t is not None and my_t>t:
                break
            if my_t<t_start:
                continue
            if i<i_start:
                continue
            if debug:
                print('# anim frame ',i,t)
            delta_t = 1/fps
            frame = []
            camera_d = None
            ichar=-1
            for character in self._characters:
                ichar+=1
                print('# get_frame()')
                f = character.get_frame(i, camera)
                print('# get_frame done()')
                if character==camera:
                    camera_d = f
                else:
                    frame.append(f)

            print('# render()')
            camera.render(frame, i, my_t, camera_d)
        camera.render_end()

@numba.njit
def interp_0_1(x0p,  x1, x1p,  t):
    D = 0
    C = x0p
    B = 3.*x1 - x1p - 2.*C - 3*D
    A = x1 - B - C - D
    return A*t*t*t + B*t*t  + C*t + D


@numba.njit
def interp(t0, x0, x0p, t1, x1, x1p, t):
    delta_x = x1 - x0
    delta_t = (t1 - t0)
    f=5.e4
    f1=f#2. #f
    f2=f#1. #f
    return x0 + interp_0_1(x0p*delta_t*f1, x1-x0, x1p*delta_t*f2,  (t - t0)/delta_t)




@numba.njit
def interp_values_chunk( t1, ids1, pos1, vel1, data1,  t2, ids2, pos2, vel2, data2, t3, ids3, pos3, data3, Ndata, i_chunk, chunks):
        N3 = len(ids3)
        ids1_count = 0
        ids2_count = 0
        for i in range(N3):
            if (i%chunks)!=i_chunk:
                continue
                
            id3 = ids3[i]
            while ids1_count<len(ids1)-1 and ids1[ids1_count]<id3:
                ids1_count +=1
            while ids2_count<len(ids2)-1 and ids2[ids2_count]<id3:
                ids2_count +=1
            id1 = ids1[ids1_count]
            id2 = ids2[ids2_count]
            if id1 == id3 and id2==id3:
                for j in range(3):
                    pos3[i][j] = 0
                    pos3[i][j] = interp(t1, pos1[ids1_count][j], vel1[ids1_count][j], t2, pos2[ids2_count][j], vel2[ids2_count][j], t3)
                for j in range(Ndata):
                    data3[i][j*2] =  data1[ids1_count][j] 
                    data3[i][j*2+1] =  data2[ids2_count][j] 
                    #pos3[i][j] =  pos1[ids1_count][j] + (t3-t1)/(t2-t1)*( pos2[ids2_count][j] - pos1[ids1_count][j] ) #, vel1[ids1_count][j], t2, pos2[ids2_count][j], vel2[ids2_count][j], t3)
            elif id1 == id3:
                for j in range(3):
                    pos3[i][j] = pos1[ids1_count,j] + (t3-t1)*vel1[ids1_count,j]*5.e4
                for j in range(Ndata):
                     data3[i][j*2] =  data1[ids1_count][j] 
                     data3[i][j*2+1] =  0.
            elif id2 == id3:
                for j in range(3):
                    pos3[i][j] = pos2[ids2_count,j] + (t3-t2)*vel2[ids2_count,j]*5.e4
                for j in range(Ndata):
                    data3[i][j*2] = 0.
                    data3[i][j*2+1] =  data2[ids2_count][j] 

@numba.njit(parallel=True)
def interp_values( t1, ids1, pos1, vel1, data1,  t2, ids2, pos2, vel2, data2, t3, ids3, pos3, data3, Ndata):
    chunks = 8
    for i_chunk in numba.prange(chunks):
        print(i_chunk)
        interp_values_chunk( t1, ids1, pos1, vel1, data1,  t2, ids2, pos2, vel2, data2, t3, ids3, pos3, data3, Ndata, i_chunk, chunks)
    


class CharacterGadget(Character):
    def __init__(self, time_factor=10., props=None):
        if props is None:
            props = []
        self.time_factor = time_factor
        super(CharacterGadget, self).__init__(props=props)
    def get_time(self, f=None, time=None):
        if f is not None:
            
            return f.header.time * self.time_factor
        else:
            return time * self.time_factor

    def get_scale_factor(self, t):

        return t/ self.time_factor
    def t_to_time(self, a):
        "formula to convert scale factor in Gyr stolen from the cosmological calculator https://home.fnal.gov/~gnedin/cc/"
        x = (1.00001-self.om0)/self.om0*a*a*a
        w = 651.3333/self.h0/np.sqrt(1.00001-self.om0)*self.q0*np.log(np.sqrt(x)+np.sqrt(x+1));
        return w*self.ureg.yr


class CameraGadgetParticlesStill(Camera, CharacterGadget):
    def __init__(self, props = None, file_format='frame%04d.png'):
        if props is None:
            props = ["GPOS","SIZE"]
        super(CameraGadgetParticlesStill, self).__init__(props=props)
        self.file_format = file_format

    def render(self, frames, i, t, camera_d):
        print('# render() ',i,t)
        f,a  = plt.subplots(1)
        gpos = camera_d['GPOS']
        size = camera_d['SIZE'].to('glength').magnitude*2
        Nfilter=100
        onion = False
        img_size = 64
        for frame in frames:
            #tit = 'KF: x%s\nPOS: %s\nID: %s'%(frame['is_keyframe'] if 'is_keyframe' in frame else '',str(frame['POS '][0]), str(frame['ID  '][0]))e
            time_a = self.get_scale_factor(frame['t'])
            time_gyr = self.t_to_time(time_a).to('yr').magnitude
            tit = 'a = %.3f, z=%.1f, %s'%(time_a, 1/time_a-1.,  'kf'  if 'is_keyframe' in frame else ' ')
            print(tit)
            a.set_title(tit)

            if not self.dmo:
                final_image_t = g3maps.smac_gen_image_matrix(img_size)
                PROPERTY = 'TEMP'
                frame_mask=frame['PTYPE']==0
                T = frame["TEMP"][frame_mask]
                T[T>1e6]=0.
                T[T<0]=0.
                if not self.dmo:
                    d = {'MASS':frame['MASS'][frame_mask], "POS ":frame["POS "][frame_mask], "RHO ":frame["RHO "][frame_mask]*self.ureg.gmass/(self.ureg.clength**3), "TEMP":T*self.ureg.dimensionless}
                else:
                    d = {'MASS':frame['MASS'][frame_mask], "POS ":frame["POS "][frame_mask]}
                
                print(d)
                g3maps.smac_add_chunk(d, 0, gpos, 2, camera_d['SIZE'], camera_d['SIZE'], img_size, PROPERTY, self.ureg.dimensionless, final_image_t, sph=False, image_w = True)
                final_image = np.nan_to_num(final_image_t)
                print(final_image)
                final_image[final_image>1e30]=0.
                final_image[final_image<0]=0.
                print('goin to imshow')
                
                we = a.imshow(final_image, extent = [gpos[0].to('glength').magnitude-size, gpos[0].to('glength').magnitude+size, gpos[1].to('glength').magnitude-size, gpos[1].to('glength').magnitude+size], vmin=0., vmax=3.e5)
                f.colorbar(we)
                print('imshow done')

                #continue
            if self.dmo:
                frame_mask=(frame['ID  ']%Nfilter==0)
                a.scatter(frame['POS '][:,0].to('glength').magnitude[frame_mask], frame['POS '][:,1].to('glength').magnitude[frame_mask],marker=',',lw=0,s =1,alpha=.1, color='black', label='now')

            else:
                frame_mask=(frame['ID  ']%Nfilter==0) & (frame['PTYPE']==4)
                a.scatter(frame['POS '][:,0].to('glength').magnitude[frame_mask], frame['POS '][:,1].to('glength').magnitude[frame_mask],marker=',',lw=0,s =1,alpha=.5, color='lightblue', label='now')
            print('# scattero')
            if onion:
                if('kf1' in frame and frame['kf1']!=None):
                    frame_mask=frame['kf1']['ID  ']%Nfilter==0
                    a.scatter(frame['kf1']['POS '][:,0].to('glength').magnitude[frame_mask], frame['kf1']['POS '][:,1].to('glength').magnitude[frame_mask], marker=',',lw=0,s =1,alpha=0.1,color='blue',label='past')
                
                if('kf2' in frame and frame['kf2']!=None):
                    frame_mask=frame['kf2']['ID  ']%Nfilter==0
                    a.scatter(frame['kf2']['POS '][:,0].to('glength').magnitude[frame_mask], frame['kf2']['POS '][:,1].to('glength').magnitude[frame_mask], marker=',',lw=0,s =1,alpha=0.1,color='red', label='future')

            a.set_xlim([gpos[0].to('glength').magnitude-size, gpos[0].to('glength').magnitude+size])
            a.set_ylim([gpos[1].to('glength').magnitude-size, gpos[1].to('glength').magnitude+size])
            #a.legend()
        print('# savefig()')
        f.savefig(self.file_format%(i))
        print('# done()')

@g3matcha.memoize
def memo_read_real(path, pos, size, blocks, ptypes, use_cache = False, dmo=False, sort=True):
    data= g3.read_particles_in_box(path, pos, size, blocks, ptypes)
    if sort:
        ids1 = data['ID  ']
    
        N = len(ids1)
    if not dmo:
        Ngas = len(data['RHO '])
        avg_rho = np.mean(data['RHO '])

        data['RHO '] = np.concatenate((data['RHO '], np.zeros(N-Ngas)+avg_rho))
        data['TEMP'] = np.concatenate((data['TEMP'], np.zeros(N-Ngas)+avg_rho))

    if sort:
        print('# begin argsort')
        ids1_argsort = np.argsort(ids1)
        for k in data:
            data[k] = data[k][ids1_argsort]
            
        print('# begin filter')

        data['DIST'] = g3.to_spherical(data['POS '], pos).T[0]

    return data

class CameraGadgetParticlesFollow(CameraGadgetParticlesStill):
    def __init__(self, units = None, ureg= None, file_format = 'frame%04d.png', dmo=False):
        super(CameraGadgetParticlesFollow, self).__init__(props = ["GPOS","SIZE","GVEL"], file_format = file_format)
        self.om0 = .272
        self.dmo = dmo
        self.q0 = self.om0 /2.
        self.h0 = .704 #TODO: read hubble constant h0 from simulation data
        self.ureg=  ureg
        self.units = units
        self.cache = {}
        #self.file_format= file_format
        
    def interp(self, kf1,kf2,t):
        pos_u =  kf1['GPOS'].to('kpc').units
        a1=self.get_scale_factor(kf1['t'])
        a2=self.get_scale_factor(kf2['t'])
        pos1 = kf1['GPOS'].to('kpc').magnitude
        pos2 = kf2['GPOS'].to('kpc').magnitude
        vel_u = kf1['GVEL'].to('kpc/s').units
        vel1 = kf1['GVEL'].to('kpc/s').magnitude/np.sqrt(a1)
        vel2 = kf2['GVEL'].to('kpc/s').magnitude/np.sqrt(a2)
        t1 = self.t_to_time(a1).to('s').magnitude
        t2 = self.t_to_time(a2).to('s').magnitude
        t3  = self.t_to_time(self.get_scale_factor(t)).to('s').magnitude        
        pos3 =[interp(t1, pos1[i], vel1[i], t2, pos2[i], vel2[i], t3) for i in range(3)]
        #pos3 = pos1 + (t3-t1)*(pos2-pos1)/(t2-t1)
        vel3  =  vel1 + (t3-t1)*(vel2-vel1)/(t2-t1)
        size3 = kf1["SIZE"] + (t3-t1)*(kf2["SIZE"]-kf1["SIZE"])/(t2-t1)
        d = {"GPOS": np.array(pos3)*pos_u,"GVEL":vel3,"SIZE":size3}


        
        return d


@g3matcha.memoize
def get_snap_info(snap_path, use_cache = False):
        f = g3.GadgetFile(g3.get_one_file(snap_path))
        return f.header.HubbleParam, f.header.time
class CharacterGadgetParticles(CharacterGadget):
    def __init__(self, blocks=None, time_factor=10., ureg=None, units=None, dmo=False):
        blocks_default = ["POS ","VEL ","MASS","ID  ","RHO ",'TEMP']
        if dmo:
            blocks_default = ["POS ","VEL ","MASS","ID  "]

        super(CharacterGadgetParticles, self).__init__(time_factor=time_factor)
        self.blocks = list(set(blocks + blocks_default)) if blocks is not None else blocks_default #we get uniq values
        self.prop = "MASS"
        self.dmo = dmo
        self.om0 = .2
        self.q0 = .1
        self.h0 = .7 #TODO: read hubble constant h0 from simulation data
        self.ureg=  ureg
        self.units = units
        self.f = None
        self.factor = None
        self.read_factor = 8.
    def add_snapshot(self, snap_path):
        h,t = get_snap_info(snap_path, use_cache = True)

        if self.ureg is None : #we still ahvent set units
            self.units = g3u.Units();
            self.units.set_context(h, 1., debug=True) #we will fo stuff in comoving so a=1.
        return self.add_keyframe(self.get_time(time=t), {"path":snap_path})
        
    # we cache reading of data
    def memo_read(self,path, pos, size, blocks, ptypes, use_cache = False):
        if self.f is None:
            self.f = g3.GadgetFile(g3.get_one_file(path))
        if self.factor is None:
            self.factor = g3u.gen_factor(self.units, self.f, blocks)
            self.factor['DIST'] = self.ureg.glength
            if not self.dmo:
                self.factor['RHO '] = self.ureg.dimensionless
        print('# memo_read_real')
        d = memo_read_real(path, pos.to('glength').magnitude, size.to('glength').magnitude, blocks, ptypes, use_cache = True, dmo = self.dmo)
        v = dict(d)
        print('#add units')
        vu =  g3u.add_units_blocks(v, blocks, self.factor, _Q=self.ureg.Quantity)
        print('# done')
        return vu
        
    def default_dict_to_frame(self, d, camera):
        print('# dict to frame')
        super(CharacterGadgetParticles, self).default_dict_to_frame(d,camera)
        if self.dmo:
            print('# dmo read')

            r = self.memo_read(d['path'], camera["GPOS"], camera["SIZE"]*self.read_factor, self.blocks, -1, use_cache=  True)
            print('# dmo readed')

        else:
            r = self.memo_read(d['path'], camera["GPOS"], camera["SIZE"]*self.read_factor, self.blocks, [0,4], use_cache=  True)

        print('# copy')
        data = dict(r)


        for k in data:
            data[k] = data[k].copy()
        size = camera["SIZE"]
        #x_distance_mask =  data['DIST']< size*self.read_factor
        #for k in data:
        #    data[k] = data[k][x_distance_mask]
        data['t'] = d['t']
        data['i'] = d['i']
        print('# done')


        return data
    

    def interp(self, kf1,kf2,t3):
        print('# begin interp')
        t1 = 0
        t2 = 0
        t3_old = t3
        ids1 = np.array([])
        ids2 = np.array([])
        pos1 = np.array([])
        pos2 = np.array([])
        vel1 = np.array([])
        vel2 = np.array([])
        if kf1 is not None:
            a1=self.get_scale_factor(kf1['t'])
            t1 = self.t_to_time(a1).to('s').magnitude
            ids1 = kf1['ID  ']
            pos1 = kf1['POS '].to('kpc').magnitude
            vel1 = kf1['VEL '].to('kpc/s').magnitude/np.sqrt(a1)
        if kf2 is not None:
            a2=self.get_scale_factor(kf2['t'])
            t2 = self.t_to_time(a2).to('s').magnitude
            ids2 = kf2['ID  ']
            pos2 = kf2['POS '].to('kpc').magnitude
            vel2 = kf2['VEL '].to('kpc/s').magnitude/np.sqrt(a2)
        from scipy import optimize
        t3  = self.t_to_time(self.get_scale_factor(t3)).to('s').magnitude
        #t3 = optimize.root(lambda t: self.t_to_time(self.get_scalefactor(t3)).to('s').magnitude - t, t1 if t1>0 else t2, method='hybr').x[0]

        l1 = len(kf1['PTYPE'])
        data1 = np.zeros((l1,4))
        data1[:,0] = kf1['PTYPE']
        if not self.dmo:
            data1[:,1] = kf1['RHO ']
            data1[:,2] = kf1['TEMP']
        data1[:,3] = kf1['MASS']

        l2 = len(kf2['PTYPE'])
        data2 = np.zeros((l2,4))
        data2[:,0] = kf2['PTYPE']
        if not self.dmo:
            data2[:,1] = kf2['RHO ']
            data2[:,2] = kf2['TEMP']
        data2[:,3] = kf2['MASS']
        rho3  = None

        if len(ids1)>0 and len(ids2)>0:
            assert(t1<t2)
            assert(t3<t2)
            assert(t1<t3)
            ids3 = np.sort(np.unique(np.sort(np.concatenate((ids1, ids2)))))
            N3 = len(ids3)
            pos3 = np.zeros(N3*3).reshape((N3,3))
            data3 = np.zeros((N3,4*2))
            print('# call numba interp')

            interp_values(t1, ids1, pos1, vel1, data1, t2, ids2, pos2, vel2, data2, t3, ids3, pos3, data3, 4)
            print('# called numba interp')

            ptype3 = np.array(data3[:,0*2], dtype=np.int32)
            rho3 = data3[:,1*2]
            temp3 = data3[:,2*2]+(data3[:,2*2+1]-data3[:,2*2])*(t3-t1)/(t2-t1)
            mass3 = data3[:,3*2]
            mass3[mass3==0.] = data3[:,3*2+1][mass3==0.]
        elif len(ids1)>0:
            pos3 = pos1 +  (t3-t1)*vel1
            ids3 = ids1
            ptype3 = np.array(data1[:,0], dtype=np.int32)
            rho3 = data1[:,1]
            temp3 = data1[:,2]
            mass3 = data1[:,3]
        elif len(ids2)>0:
            pos3 = pos2 +  (t3-t2)*vel2
            ids3 = ids2
            ptype3 = np.array(data2[:,0], dtype=np.int32)
            rho3 = data2[:,1]
            temp3 = data2[:,2]
            mass3 = data2[:,3]

        else:
            print('kf1')
            print(kf1)
            print('kf2')
            print(kf2)
            print(t3)
            raise Exception('?')
        pos3 = self.ureg.Quantity(pos3, self.ureg.kpc)
        print('# interped')

        lego= {"POS ": pos3, "ID  ":ids3,"kf1":kf1, "kf2":kf2, "t":t3_old,"PTYPE":ptype3, 'RHO ':rho3, 'TEMP':temp3, 'MASS':self.ureg.Quantity(mass3, self.ureg.gmass)}

        print('# goin to return')
        return lego
