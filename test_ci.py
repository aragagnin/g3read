import g3read, g3matcha, numpy as np, urllib.request, os



def download_from_googledrive(url_googledrive, outname):
    file_id = url_googledrive.split('/')[-2]
    download_url =  'https://docs.google.com/uc?export=download&id=' + file_id
    print('downloading...', download_url,'into',outname)
    urllib.request.urlretrieve(download_url, outname)

if not os.path.isfile('snap_024'):
    download_from_googledrive('https://drive.google.com/file/d/1b4u1HLaXyBIXLfnqOYhZsA_iAWme3KWA/view?usp=sharing', 'snap_024')

array = np.array
uint32 = np.uint32
int32 = np.int32
float32 = np.float32

#
# test 1: read of header
#
expected_header = {'mass': array([0.        , 0.26337509, 1.3076342 , 0.        , 0.        ,
       0.        ]), 'time': 0.6801823040103018, 'redshift': 0.4701940848858872, 'flag_sfr': 0, 'flag_feedback': 0, 'flag_cooling': 0, 'num_files': 1, 'BoxSize': 48000.0, 'Omega0': 0.272, 'OmegaLambda': 0.728, 'HubbleParam': 0.704, 'flag_stellarage': 0, 'flag_metals': 0, 'flag_entropy_instead_u': 0, 'flag_doubleprecision': 0, 'flag_ic_info': 3, 'lpt_scalingfactor': 0.0, 'endian': '=', 'npart': array([     0, 531441, 531441,      0,      0,      0], dtype=uint32), 'npartTotal': array([     0, 531441, 531441,      0,      0,      0], dtype=int32), 'NallHW': array([0, 0, 0, 0, 0, 0], dtype=int32)}
f = g3read.GadgetFile('snap_024')

print(f.header.__dict__) # use this print to obtain the `expected_header`

for k,v in expected_header.items():
    hv = f.header.__dict__[k]
    print(k,':',v,'=',hv)
    assert(np.all((hv-v)<1e-7) if  isinstance(hv, np.ndarray) else hv==v)


#
# test 2: read of a block
#

vel0 = array(
    [-170.54239,   -67.39497,   -46.011555]
    , dtype=np.float32)
velm1 = array(
     [-185.24017,  -79.13706,   -60.12431 ]
    , dtype=np.float32)
vel = g3read.read_new('snap_024', 'VEL ', 1)
print(vel0)
print(vel)
print(velm1)
for i in range(3):
    err = vel[0][i]-vel0[i]
    if(np.abs(err)>1e-10):
        raise Exception("wrong velocity read error:%f"%err)
    err = vel[-1][i]-velm1[i]
    if(np.abs(err)>1e-10):
        raise Exception("wrong velocity read error:%f"%err)



#
# test 2: read of a block OOP way
#
f = g3read.GadgetFile('snap_024')
vel = f.read_new('VEL ',1)
print(vel0)
print(vel)
print(velm1)
for i in range(3):
    err = vel[0][i]-vel0[i]
    if(np.abs(err)>1e-10):
        raise Exception("wrong velocity read error:%f"%err)
    err = vel[-1][i]-velm1[i]
    if(np.abs(err)>1e-10):
        raise Exception("wrong velocity read error:%f"%err)

   
#
# test 3: read of subfind haloes
#
first_halo = {'GPOS': array([ 6322.5938, -4340.9375, -2297.6562], dtype=float32), 'RCRI': np.float32(1551.3738), 'MCRI': np.float32(41838.207), 'GOFF': np.int32(0), 'GLEN': np.int32(1118201), 'ihalo': 0, 'ihalo_in_file': 0, 'ihalo_file': 'sub_060.0', 'boxsize': 1000000.0}
glen0 = first_halo['GLEN']
goff0 = first_halo['GOFF']
for halo in g3matcha.yield_haloes('sub_060.0'):
    print(halo) # use it to produce `first_halo`
    for k,v in first_halo.items():
        hv = halo[k]
        assert(np.all(hv==v) if  isinstance(hv, np.ndarray) else hv==v)

    break


#
# test 3: read of subfind subhaloes
#
first_subhalo = {'SMST': array([4.1611465e+03, 3.3093648e+04, 0.0000000e+00, 0.0000000e+00,
       7.1197760e+02, 1.6695596e+15], dtype=float32), 'SPOS': array([ 6322.5938, -4340.9375, -2297.6562], dtype=float32), 'SOFF': np.int32(0), 'SLEN': np.int32(879366), 'GRNR': 0, 'ihalo': 0, 'ifile': 0, 'isubhalo': 0}
for subhalo in g3matcha.yield_subhaloes('sub_060.0',0):
    print(subhalo) # use it to produce `first_halo`
    for k,v in first_subhalo.items():
        hv = subhalo[k]
        assert(np.all(hv==v) if  isinstance(hv, np.ndarray) else hv==v)

    break


#
# test 4 get halo IDs
#
ids, ifile, goff_file = g3matcha.get_halo_ids('sub_060.0', goff0, glen0, ifile_start=0, goff_start=0)
ids.dtype = np.uint32
assert(ids[1] == 2149688460)
assert(ifile == 0)
assert(goff_file == 959471)

print('\n\nAll test passed!\n\n')
