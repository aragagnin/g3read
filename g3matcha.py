#!/usr/bin/env python
"""
These routines help match Gadget haloes between boxes and snapshots; to loop over haloes and subhaloes of a simulation and to match/extract subfind IDs.

Antonio Ragagnin (2020) <ragagnin.antonio@gmail.com>

"""
try:
        from . import g3read as g3
except:
        import g3read as g3
        
import numpy as np, yaml, json, sys
from collections import OrderedDict

import os
import copy

import pickle; #_pickle as cPickle



class pdict(OrderedDict): #persistent dictionary
        def __init__(self, filename = None, *args, **kw):
            super(pdict,self).__init__(*args, **kw)
            self.filename = filename
            if filename is not None:
                if os.path.isfile(filename):
                    self.restore()
        def restore(self):
            if self.filename is None:
                raise Exception('cannot restore if filename is None')
            if cache_type == 'pickle':
                if os.path.isfile(self.filename):
                        with open(self.filename,'rb') as f:
                                d = pickle.load(f)
                        self.update(d)
                 
        
        def store(self):
            if self.filename is None:
                raise Exception('cannot store if filename is None')

            elif cache_type == 'pickle':
                with open(self.filename+'~','wb') as f:
                        pickle.dump(dict(self),f)
                os.rename(self.filename+'~',self.filename)

#
# here below a small cache system to parse many FoF files fast.
#

class LimitedSizeDict(pdict):
    #stolen from https://stackoverflow.com/questions/2437617/how-to-limit-the-size-of-a-dictionary
    def __init__(self, filename = None, size_limit=20):
        self.size_limit = size_limit
        super(LimitedSizeDict, self).__init__(filename = filename)

    def trim(self):
        while len(self) > self.size_limit:
            self.popitem(last=False)

def dict_to_pairs(d):
    return [(i, d[i]) for i in d]

debug = False
cache_from_filename_only = True
cache_filename = None
cache = None
size_limit = 2000
cache_type = 'pickle'

def memoize(func):
    global cache, size_limit
    """
    decorator to cache function results, taken from: https://dbader.org/blog/python-memoization
    """
    def memoized_func(*args, **kw):
        global cache, size_limit
        if 'use_cache' not in kw or kw['use_cache']==False:
                return func(*args, **kw)
        else:
            if cache is None:
                cache = LimitedSizeDict(filename = cache_filename, size_limit = size_limit)
                if cache_filename is not None:
                        if debug:
                                print('# prima di restore ', cache.keys())
                        cache.restore()
                        if debug:
                                print('# dopo di restore ', cache.keys())
            funcname = func.__name__ if cache_from_filename_only else str(func)
            k = str((funcname, tuple(args), tuple(dict_to_pairs(kw))))
            if k in cache.keys():
                return copy.deepcopy(cache[k])
            if debug:
                print('->',k)
            result = func(*args, **kw)
            cache[k] = result

            cache.trim()
            if cache_filename is not None:
                    cache.store()

            return result

    return memoized_func

#
# now functions with @memoize will cache results
#


@memoize
def get_fof_file(filename, use_cache=False):
    return g3.GadgetFile(filename, is_snap=False)

@memoize
def read_new(filename, block, ptypes, use_cache=False):
    s = get_fof_file(filename, use_cache = use_cache)
    return s.read_new(block, ptypes)

def numpy_to_dict(data, blocks):
    """
    super-compact creation of array of dicts of read_new data
    you can access data as res[0]['MASS'] instead of res['MASS'][0] as returned by read_new
    """
    keys = list(data.keys())
    k0 = blocks[0]
    result = []
    for i in range(len(data[k0])):
        d = {}
        for block in g3.iterate(blocks):
                d[block] = data[block][i]
        result.append(d)
    return result

@memoize
def read_new_dict(filename, blocks, ptypes, use_cache = False, _filter = None):
        read_new_result = read_new(filename, blocks, ptypes, use_cache = use_cache)
        
        result =  numpy_to_dict(read_new_result, blocks)
        keys =  [k for k in read_new_result.keys()]
        for key in keys:
                del read_new_result[key]
        return result

def get_halo_ids(groupbase, goff, glen, ifile_start=0, goff_start=0, use_cache = False):
    """
    given the basepath of FoF/SubFind output in `groupbase`, a goff and glen, returns an array of IDs, the last file where it found it and the offset of said file.
    The last two numbers are very useful when saving IO time in reading ids of increasing haloes ids. 

    For instance
    
        halo0_ids, ifile_start, goff_start =  get_halo_ids(groupdase, halo0_goff, halo0_glen)
        halo1_ids, ifile_start, goff_start =  get_halo_ids(groupdase, halo1_goff, halo1_glen, ifile_start = ifile_start, goff_start = goff_start)
        halo2_ids, ifile_start, goff_start =  get_halo_ids(groupdase, halo2_goff, halo2_glen, ifile_start = ifile_start, goff_start = goff_start)

    is much faster than the following which will start reading files from beginning at every call:

        halo0_ids =  get_halo_ids(groupdase, halo0_goff, halo0_glen)[0]
        halo1_ids =  get_halo_ids(groupdase, halo1_goff, halo1_glen)[0]
        halo2_ids =  get_halo_ids(groupdase, halo2_goff, halo2_glen)[0]
    
    """
    
    ifile=-1
    partial_ids = None
    finish=False
    goff_file = goff_start

    if debug>0:
        print('# get_halo_ids: read files, ifile_start= ',ifile_start) 
    for group_file in g3.yield_all_files(groupbase):
        ifile+=1
        if debug>0:
            print('# get_halo_ids: iterating file ',group_file)
        if(ifile < ifile_start):
            continue
        if debug>0:
            print('# get_halo_ids: read file ',group_file)
        ids_in_file = read_new(group_file,  'PID ', 2, use_cache = use_cache)
        glen_file = len(ids_in_file)
        if debug>0:
            print( '# get_halo_ids:' ,group_file, 'halo goff',goff, 'halo glen',glen, 'file cumulative goff', goff_file,' file glen',glen_file)
        if True:
            if (goff+glen<=goff_file):
                if debug>0:
                    print( '# get_halo_ids:   (goff+glen<=goff_file) => we have read all IDs')
                finish = True
            elif goff>(goff_file+glen_file): 
                if debug>0:
                    print( '# get_halo_ids:  goff>(goff_file+glen_file) => we didnt reach our first file yet' )
               
            else:
                # papabile
                start_reading = goff - goff_file
                if start_reading<0:
                    if debug>0:
                            print( '# get_halo_ids: IDs started in the prev file, I pick IDs from the begginning of the file')
                    start_reading = 0
                
                #end_reading = start_reading+glen
                end_reading = goff - goff_file + glen # CRF 21/04/2021
                if end_reading>glen_file:
                    if debug>0:
                            print( '# get_halo_ids: IDs will finish in the next file(s): I read up to end of file')
                    end_reading = glen_file
                else:
                    if debug>0:
                            print( '# get_halo_ids: Since we do not read till the end of file => this is the last file to read')
                    finish = True
                if debug>0:
                        print( '# get_halo_ids: I read in the following range, ', start_reading, end_reading)
                

                _partial_ids = ids_in_file[start_reading:end_reading]
                del ids_in_file

                
                if partial_ids is None:
                    partial_ids = _partial_ids
                else:
                    partial_ids = np.concatenate((partial_ids, _partial_ids))
            if debug>0:
                print('# get_halo_ids: partial_ids.shape = ',partial_ids.shape)    

        if debug>0:
                print( '#')


        if finish :
            break
        goff_file += glen_file        
    return (partial_ids, ifile, goff_file)

def yield_haloes(groupbase, ihalo_start=0, ihalo_end=None, min_value=None, min_block_value =None, use_cache = False, blocks=None, with_ids=False):
    """
    returns all haloes (each in a dict) and their FoF data (GLEN, GPOS, RCRI, MCRI) 
    given the path of catalog groupbase.
    
    You can set when to stop the research by givin the minimum value of a given block. 
    E.g. min_value=10000, min_block_value='GLEN' will stop after finding a halo with less than min_value GLEN.
    """
    icluster = -1
    ifile1 = -1
    
    if blocks is None:
        blocks = ('GPOS' , 'RCRI', 'MCRI')
    if 'MCRI' not in blocks:  blocks = blocks+('MCRI',)
    if 'GOFF' not in blocks:  blocks = blocks+('GOFF',)
    if 'GLEN' not in blocks:  blocks = blocks+('GLEN',)
    if min_block_value is not None and min_block_value not in blocks:  blocks = blocks+(min_block_value,)
    if with_ids:
        ifile_ids_start = 0
        ifile_ids_goff = 0
    for group_file in g3.yield_all_files(groupbase):

        s = get_fof_file(group_file, use_cache = use_cache)
        nclusters_in_file = s.header.npart[0]
        clusters_in_file = read_new_dict(group_file,  blocks, 0, use_cache = use_cache)
        boxsize1 = s.header.BoxSize
        if min_value is not None  and   (clusters_in_file[0][min_block_value] < min_value):
            return
        for icluster_file in range(nclusters_in_file):

            icluster = icluster+1
            
            if icluster<ihalo_start:
                continue
            if ihalo_end is not None and icluster>ihalo_end:
                return
            halo_info = dict(clusters_in_file[icluster_file])
            if with_ids:
                goff = halo_info['GOFF']
                glen = halo_info['GLEN']
                halo_ids, ifile_ids_start, ifile_ids_goff =  get_halo_ids(groupbase, goff, glen, ifile_start = ifile_ids_start, goff_start = ifile_ids_goff, use_cache = use_cache)
                halo_info['ids'] = halo_ids
                
            
            halo_info['ihalo'] = icluster
            halo_info['ihalo_in_file'] = icluster_file
            halo_info['ihalo_file'] = group_file
            #this prop  is not really related to the halo but I need it for periodic distance
            halo_info['boxsize'] = boxsize1 
            
            yield halo_info
            del halo_info
            
@memoize
def get_all_files(groupbase, use_cache=False):
    return list(g3.yield_all_files(groupbase))

    
def yield_subhaloes(groupbase, ihalo, ifile_start=None,  use_cache = False, blocks=None, with_ids = False, halo_ids=None, halo_goff=0):
    """
    returns all subhaloes (each in a dict) and their FoF data (GLEN, GPOS, RCRI, MCRI) 
    given the path of catalog groupbase
    """
    icluster = -1
    ifile1 = -1
    if blocks is None:
        blocks = ('SMST','SPOS','SOFF')
    if 'SOFF' not in blocks:   blocks = blocks+('SOFF',)
    if 'SLEN' not in blocks:   blocks = blocks+('SLEN',)

    found_first_subhalo = False
    isubhalo = -1
    ifile=-1

    for group_file in  get_all_files(groupbase, use_cache=use_cache):
        ifile+=1

        if ifile_start!=None and ifile<ifile_start:
            continue
        grnr = read_new(group_file, 'GRNR', 1, use_cache = use_cache)
        #minimum halo associated to subhaloes in this file is greater than what we seek
        if np.min(grnr)>ihalo:
            return
        data = read_new(group_file, blocks, 1, use_cache = use_cache)
        
        if not np.any(grnr==ihalo):
            
            #this file dosent contain subhaloes of this halo
            if found_first_subhalo: 
                #we read all files with this subhalo data
                return
            continue


        found_first_subhalo = True
        subhaloes_in_file = numpy_to_dict(data, blocks)
        #print('loppi')
        for isubhalo_in_file, subhalo in enumerate(subhaloes_in_file):
            if (grnr[isubhalo_in_file]!=ihalo):
                continue
            isubhalo+=1
            subhalo['GRNR'] = grnr[ihalo]
            subhalo['ihalo'] = ihalo
            subhalo['ifile'] = ifile
            subhalo['isubhalo'] = isubhalo
            if with_ids:
                subhalo['ids'] = halo_ids[subhalo['SOFF']-halo_goff:subhalo['SOFF']+subhalo['SLEN']-halo_goff]
            yield subhalo


def yield_matches(snapbase1, gpos1, r200c1, groupbase2, snapbase2,  ids_block1, ids_block2, min_value, min_ids_len_factor, max_distance, max_r200c_factor_ids, min_block_value='MCRI', use_cache = False):
    """
    loop over haloes and find matches with the halo on gpos1 and radius r200c1 in simulation snapbase1.
    In simulation1 read block ids_block1, in simulation2 read block ids_block2 - bc DMO runs  have blocks shifted wrt the BAO counterpart.

    """
    ids1 = None #cache ids list of this halo
    icluster2 = -1

    for cluster2 in yield_haloes(groupbase2, 0, min_value=min_value, min_block_value=min_block_value, use_cache = use_cache):
        glen2  = cluster2['GLEN']
        gpos2  = cluster2['GPOS']
        r200c2  = cluster2['RCRI']
        m200c2  = cluster2['MCRI']
        boxsize1  = cluster2['boxsize']

        distance = g3.periodic_distance(gpos1, gpos2, periodic=boxsize1)
        cluster2['distance'] = distance
        cluster2['int_frac'] = np.nan


        if(distance<=max_distance):

            if (ids_block1 is not None) and (ids_block2 is not None):
                if ids1 is None: #check if halo ids are in cache or if we must obtain them

                    dm_data1  = g3.read_particles_in_box(snapbase1, gpos1,  r200c1, ['POS ', 'ID  '],  ids_block1)
                    mask1 = g3.to_spherical(dm_data1['POS '], gpos1).T[0] < (r200c1 * max_r200c_factor_ids)
                    ids1 = np.sort(dm_data1['ID  '][mask1])
                    
                dm_data2  = g3.read_particles_in_box(snapbase2, gpos2,  r200c2, ['POS ', 'ID  '],  ids_block2)
                mask2  = g3.to_spherical(dm_data2['POS '], gpos2).T[0] < (r200c2 * max_r200c_factor_ids)
                ids2 =  np.sort(dm_data2['ID  '][mask2])
                ids_int  = np.intersect1d(ids1, ids2)
                int_len = len(ids_int)


                int_frac = float(int_len)/float(len(ids1))
                cluster2['int_frac'] = int_frac

                if int_frac > min_ids_len_factor:

                    yield cluster2
                    return
            else:
                yield cluster2
                return
                


def nfw_fit_fast_cu(mass,rs,R,nbins=50):
    import scipy
    import scipy.optimize

    r_bins = np.logspace(np.log10(R/nbins),np.log10(R), nbins)
    mass_m,mass_bin = np.histogram(rs, bins=r_bins, weights=mass)
    r_r,r_bin = np.histogram(rs, bins=r_bins, weights=rs)
    r_n,r_bin = np.histogram(rs, bins=r_bins)
    r_avg = r_r/r_n
    rho_my = mass_m/(4.*np.pi*(r_bin[1:]**3-r_bin[:-1]**3)/3.)

    distance_cu_to_comcgs = 3.086e21
    mass_cu_to_comcgs = 1.99e33*1.e10
    distance3_cu_to_comcgs = distance_cu_to_comcgs**3
    density_cu_to_comcgs = mass_cu_to_comcgs/distance3_cu_to_comcgs
    myprofilo_nfw=lambda r,rho0,rs: rho0 / ( (r/rs) * ((1.+r/rs)**2.))
    #print(rho_my*density_cu_to_comcgs*1e24, r_avg/R)
    minimize_me = lambda x: np.sqrt(

        np.sum(
            np.abs(
                np.log10(myprofilo_nfw(r_avg/R,x[0],x[1])/(rho_my*1e24*density_cu_to_comcgs))
                )**2
            ))


    x0=[0.05,0.5]
    method='L-BFGS-B'
    xbnd=[[0.001,50.0],[0.01,10.0]]
    r=scipy.optimize.minimize(minimize_me,x0,method=method,bounds=   xbnd)
    return {
        "rho0":r.x[0]/(1e24*density_cu_to_comcgs),
        "c":1./r.x[1]
    }




#
# HERE BELOW ONLY TEST ROUTINES
#

def printf(s,e=False):
    fd=sys.stderr if e else sys.stdout
    fd.write(s)




def test(snapbase1, groupbase1, ids_block1, snapbase2, groupbase2, ids_block2,spaces="", printf=printf, nhalos=10, use_cache=False) :
    """
    run matcha algorithm between snapshot snapbase1, groupbase1, ids_block1 and  snapbase2, groupbase2, ids_block2.
    for Magenticum: set ids_block = 1 for BAO sims and set ids_block=2 for DMO sims. This because Magneticum DMO have matter blocks shifted of 1.
    """
    ihalo_start1 = 0
    ihalo_end1 =  nhalos
    max_mcri_factor = .1 #match only haloes with at least .1 of mass of the original halo
    max_distance = 1.e3 #match only haloes within 1Mpc/h
    max_r200c_factor_ids = 0.1  #matches IDs within 0.1R200c 
    min_ids_len_factor = 0.3 #we want at least 30% ofthe IDs to be in the other halo


    
    #loop over first snapshot
    for cluster1 in yield_haloes(groupbase1, ihalo_start1, ihalo_end1, use_cache=use_cache):
        printf(spaces+"- ihalo: %d\n"%(cluster1['ihalo'])+
               spaces+"  glen: %d\n"% (cluster1['GLEN'])+
               spaces+"  center: [%.1f, %.1f, %.1f]\n"%( cluster1['GPOS'][0], cluster1['GPOS'][1], cluster1['GPOS'][2])+
               spaces+"  r200c: %.2f\n"%(cluster1['RCRI'])+
               spaces+"  M200c: %.2f\n"% (cluster1['MCRI'])+
               spaces+"  matches: \n")
        

        for cluster2 in yield_matches(snapbase1,  cluster1['GPOS'],  cluster1['RCRI'],
                                      groupbase2, snapbase2,   ids_block1,  ids_block2,
                                      max_mcri_factor* cluster1['MCRI'] ,  min_ids_len_factor,   max_distance, max_r200c_factor_ids,
                                      use_cache=use_cache):
            #print(cluster2)
            printf(spaces+"  - ihalo: %d\n"%(cluster2['ihalo'])+
               spaces+"    glen: %d\n"% (cluster2['GLEN'])+
               spaces+"    center: [%.1f, %.1f, %.1f]\n"%( cluster2['GPOS'][0], cluster2['GPOS'][1], cluster2['GPOS'][2])+
               spaces+"    r200c: %.2f\n"%(cluster2['RCRI'])+
               spaces+"    M200c: %.2f\n"% (cluster2['MCRI'])+
                   spaces+"    distance: %.1f\n"%(cluster2['distance'])+
                   spaces+"    int_frac: %.1f\n\n"%( cluster2['int_frac']))


def test_many(list_of_snapbase_groupbase1_ids_block , spaces="", printf=printf, nhalos=10 ,use_cache = False) :

    first_simulation = list_of_snapbase_groupbase1_ids_block[0]
    other_simulations = list_of_snapbase_groupbase1_ids_block[1:]
    for other_simulation in other_simulations:
        
        printf(spaces+'- match: ["%s","%s"]\n  haloes:\n'%(first_simulation['label'] if 'label' in first_simulation else first_simulation['snapbase'],
                                                           other_simulation['label'] if 'label' in other_simulation else other_simulation['snapbase']))

        test(first_simulation['snapbase'], first_simulation['groupbase'], first_simulation['ids_block'],
             other_simulation['snapbase'], other_simulation['groupbase'], other_simulation['ids_block'],
             spaces='  ',
             nhalos=nhalos,
             use_cache = use_cache
        )
        

def main():
    try:
        snapbase1, groupbase1, ids_block1, snapbase2, groupbase2, ids_block2 = sys.argv[1:]
    except:
        printf('\nMatcha 1.0 Antonio Ragagnin (2020) <ragagnin.antonio@gmail.com>\n\n')
        printf('usage: ./matcha.py snapbase1 groupbase1 ids_block1 snapbase2 groupbase2 idsblock2\n',e=True)
        printf('example: ./matcha.py /HydroSims/Magneticum/Box4/uhr_test/snapdir_136/snap_136 '
               '/HydroSims/Magneticum/Box4/uhr_test/groups_136/sub_136 1 '
               '/HydroSims/Magneticum/Box4/uhr_dm/snapdir_136/snap_136 '
               '/HydroSims/Magneticum/Box4/uhr_dm/groups_136/sub_136 2\n\n',e=True)
        sys.exit(1)
    test(snapbase1, groupbase1, int(ids_block1), snapbase2, groupbase2, int(ids_block2), use_cache=True)        

if __name__ == "__main__":
    main()
                
