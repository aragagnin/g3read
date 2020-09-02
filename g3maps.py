#!/usr/bin/env python3
# g3maps is inspired by Klaus' SMAC (https://wwwmpa.mpa-garching.mpg.de/~kdolag/Smac/) and it produce 
# 2D syntetic observations from Gadget snapshots. The advantages w.r.t. SMAC are that (1) each input is given  
# given with respective units, and (2) it can produce maps of objects that do not fit the RAM capacity of the machine.
# Units are handled with the package `pint` (https://pint.readthedocs.io/en/stable/).
#
# by Antonio Ragagagnin <antonio.ragagnin@inaf.it>
#

import g3read_units as g3u
import g3read as g3
import sys
import numpy as np
from numba import jit


# the input values are given with a file or a dict and
# `xyz_keys` stores the key name of the coordinates of center of the map
xyz_keys = ['CENTER_X', 'CENTER_Y', 'CENTER_Z'] 

def parse_ureg_expressions(kv, to_parse, ureg):
    """Reads input values with units (e.g. '3000 kpc') and parses it in a  `pint` value).
    Input values came from the config file and are passed in a dict `kv`.
    `to_parse` is a dictionary that stores the type of the value of a given key (e.g. the number of pixels IMG_SIZE is `int`),
    and `ureg` is a `pint` instance """  
    res = {}
    for k in to_parse:
        res[k] = kv[k].evaluate(ureg)
    return res


@jit(nopython=True)
def     add_to_grid(final_image_t, masked_x, masked_y, masked_h, masked_w, bin_min, delta_bin, n_bins):
    """This routine adds a chunk of particles with sky positions `masked_x, masked_y`, 
    smoothing length  `masked_h,` and value `masked_w` (e.g. paticle mass) to a FIT buffer `finalt_image_t`.
    Data is inserted in chunks in order to be able to process objects that do not fit into memory.
    """
    
    N = len(masked_x)

    tot_max_bin_spread = 0
    for k in range(N):
        bin_i = int((masked_x[k]-bin_min)/delta_bin)
        bin_j = int((masked_y[k]-bin_min)/delta_bin)
        max_bin_spread = int(masked_h[k]/delta_bin)
        tot_max_bin_spread+=max_bin_spread
        v = masked_w[k]
        for i in range(bin_i - max_bin_spread, bin_i + max_bin_spread + 1):
            if i>=n_bins or i<0:
                continue
            for j in range(bin_j - max_bin_spread, bin_j + max_bin_spread + 1):
                if j>=n_bins or j<0:
                    continue
                final_image_t[j][i]+=v
    if N>0:
        return tot_max_bin_spread/N
    else:
        return np.nan

def smac_for_poors(kv):
    """here we read input data `kv`, read chunks of particles and send them to `add_to_grid`"""
    #set the snapshots' scalefactor and hubble factor into units.
    units = g3u.get_units(kv['SNAP_PATH'])
    #produce a set of pint units from the snapshots data
    ureg = units.get_u()
    
    kv.update(parse_ureg_expressions(kv, kv['TOPARSE'], ureg))

    
    
    gpos_cu = np.array([kv[v].to('glength').magnitude for v in xyz_keys])*ureg.glength
    gpos = np.array([kv[v].to('glength').magnitude for v in xyz_keys])*ureg.glength

    img_xy_size = kv['IMG_XY_SIZE']
    img_z_size  = kv['IMG_Z_SIZE']
    img_size = kv['IMG_SIZE']
    max_xyz_cu = max( kv['IMG_Z_SIZE'].to('glength').magnitude, kv['IMG_XY_SIZE'].to('glength').magnitude)

    print('#')
    print('# XY SIZE = ', img_xy_size.to('glength'), ' = ', img_xy_size.to('kpc'))
    print('# 1e10 Msun =  ',ureg('1e10 Msun').to('gmass'))
    print('# Output image: 2D total mass')
    print('#')
    
    n_bins=[img_size]*2
    final_image_t = np.zeros(n_bins)
    #final_image_n = np.zeros(n_bins)

    ptypes = [int(x.strip()) for x in kv['PTYPES'].split(',')]
    #print('# reading ',kv['SNAP_PATH'],' centerd on code units ',gpos,' radius= ',max_xyz_cu, 'ptypes', ptypes)
    PROPERTY = kv['PROPERTY']
    bins = [np.linspace(-.5,.5,img_size+1), np.linspace(-.5, .5, img_size+1)]
    for _res in  g3u.yield_particles_blocks_in_box(kv['SNAP_PATH'], gpos_cu, max_xyz_cu,[PROPERTY,'POS ','RHO '], ptypes,debug=True, units=units):
        for ptype in _res:
            P = res = _res[ptype]
            N_part_type = len(P[PROPERTY])

            if N_part_type == 0:
                continue

            if ptype!=0:
                P['RHO '] = np.ones(N_part_type) * ureg.parse_expression(' 1.e-2 gmass/glength^3 ')
                
            P['HSML'] = (P['MASS']/(4./3. * np.pi * P['RHO '] ))**(1./3.)

            rel_poses = P['POS '] - gpos

            axis_z=2

            pos_x = rel_poses[:,(axis_z+1)%3]
            pos_y = rel_poses[:,(axis_z+2)%3]
            pos_z = rel_poses[:,(axis_z+0)%3]
            
            norm_x = (pos_x/img_xy_size).to('').magnitude 
            norm_y = (pos_y/img_xy_size).to('').magnitude 
            norm_z = (pos_z/img_z_size).to('').magnitude 

            mask = (norm_x<=.5)&(norm_x>=-.5)&(norm_y<=.5)&(norm_y>=-.5)&(norm_z<=.5)&(norm_z>=-.5)
            
            masked_x = norm_x[mask]
            masked_y = norm_y[mask]
            masked_h = (P['HSML'][mask]/img_xy_size).to('').magnitude
            masked_w = (P[PROPERTY][mask]).to(kv['RESULT_UNITS']).magnitude

            #final_image_t += np.histogram2d(masked_x,masked_y,weights=masked_w,bins=bins)[0].T
            avg_bin_spread = add_to_grid(final_image_t, masked_x, masked_y, masked_h, masked_w, bins[0][0], bins[0][1]-bins[0][0], img_size)
            print ('# avg bin spread ', avg_bin_spread)
    final_image = np.nan_to_num(final_image_t)

    print('#')
    # write the FIT file down
    print(final_image.shape)
    import matplotlib.pyplot as plt
    from astropy.io import fits
    outfile = kv['PREFIX_OUT']+'.fits'
    print("# writing...", outfile)
    header = fits.Header()
    header.append(('EXTENSION','IMAGE','IMAGE extension'))
    header.append(('NAXIS', 2,' number of data axes'))
    header.append(('NAXIS1', kv['IMG_SIZE'] ,' length of data axis 1'))
    header.append(('NAXIS2', kv['IMG_SIZE'],' length of data axis 1'))
    header.append(('OUTTYPE','2D Total density',' outtype'))
    header.append(('OUTUNIT',kv['RESULT_UNITS'],' untis of map'))
    header.append(('PROJECT',['x','y','z'][axis_z],' direct. of proj'))
    header.append(('REDSHIFT',repr((1./ureg.parse_expression('scalefactor')-1.).to('').magnitude),'  redshift of snapshot'))
    header.append(('PSIZEPKC',(img_xy_size/img_size).to('kpc').magnitude,' kpc/pixel size'))
    header.append(('SNAPPATH',kv['SNAP_PATH'],' sim path'))

    

    hdu = fits.HDUList([
        fits.PrimaryHDU(header=fits.Header.fromstring("""\
SIMPLE  =                    T / conforms to FITS standard
BITPIX  =                    8 / array data type
NAXIS   =                    0 / number of array dimensions
EXTEND  =                    T
        """, sep='\n')

        ),
        
        fits.ImageHDU(final_image,
                      header=header
        )
    ])
    hdu.writeto(outfile, clobber=True)
    print('#')

# factory for ureg pint objects from a string
class Uregexpr(object):
    def __init__(self, expression):
        self.expression = expression
    def evaluate(self, ureg):
        return ureg.parse_expression(self.expression)
    def __str__(self):
        return 'ureg '+self.expression

# type of input paramters, Uregexpr means the input may be  a string like '1000. kpc'
types={
    "CENTER_X": Uregexpr,
    "CENTER_Y": Uregexpr,
    "CENTER_Z": Uregexpr,
    "IMG_XY_SIZE": Uregexpr,
    "IMG_Z_SIZE": Uregexpr,
    "IMG_SIZE": int,

}    
    


def parse_config_file(filename):
    " convert a INI file to a dict where values have types from `types` "

    kv = {'TOPARSE':[]}
    with open(input_file,'r') as f:
        for line0 in f:
            line1=line0
            if '#' in line0:
                line1 = line0.split('#',1)[0]
            if '=' in line1:
                k0,v0 = line1.split('=',1)
                k1=k0.strip()
                v1=v0.strip()
                v2=v1
                if k1 in types:
                    v2=types[k1](v2)
                    if types[k1].__name__=='Uregexpr':
                        kv['TOPARSE'].append(k1)
                kv[k1]=v2
                
                print('#',k1,'=',str(v2))
    return kv

#
# if we execute the script we read an input file specified in the command line
# and pass it to smac_for_poors.
#
if __name__=='__main__':
    print('#')
    print('# g3maps.py: Gadget large maps v0.1a ')
    print('#')

    if (len(sys.argv)!=2):
        
        print("Usage: ./smac_for_poors.py input_file.inp")
        print("Got: ", sys.argv)
        sys.exit(1)
        
    input_file = sys.argv[1]
    kv = parse_config_file(input_file)

    if kv['JOB']=='2DMAP':
        smac_for_poors(kv)
    else:
        print('NOT IMPLEMENTED JOB TYPE ', kv['JOB']);
        sys.exit(1);
