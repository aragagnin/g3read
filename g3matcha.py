#!/usr/bin/env python
"""
These routines help match Gadget haloes between boxes and snapshots; to loop over haloes and subhaloes of a simulation and to match/extract subfind IDs.

Antonio Ragagnin (2020) <ragagnin.antonio@gmail.com>
"""

try:
        from . import g3read as g3
except:
        import g3read as g3
        
import numpy as np,  sys
from collections import OrderedDict

import os
import copy

import pickle; #_pickle as cPickle



class pdict(OrderedDict): #persistent dictionary
        def __init__(self, filename = None, read_only = False,*args, **kw):
            super(pdict,self).__init__(*args, **kw)
            self.filename = filename
            self.read_only = read_only
            if filename is not None:
                if os.path.isfile(filename):
                    with open(filename,'rb') as f:
                        d = pickle.load(f)
                    self.update(d)
        def store(self):
            if self.filename is None:
                raise IOException('cannot store if filename is None')
            if self.read_only:
                raise IOException('dict is read only')
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
size_limit = 2000
cache_type = 'pickle'
recache = None
cache_filename_default = 'cache'
cache_filenames = {}


def memoize(func):
    global size_limit,recache
    """
    decorator to cache function results, taken from: https://dbader.org/blog/python-memoization
    """
    def memoized_func(*args, **kw):
        global size_limit, recache
        if 'use_cache' not in kw or kw['use_cache']==False:
                return func(*args, **kw)
        else:
            if isinstance(kw['use_cache'], str):
                    cache_filename = kw['use_cache']
            else:
                    cache_filename = cache_filename_default
                    
            if cache_filename not in cache_filenames or recache==True:
                cache_filenames[cache_filename] = LimitedSizeDict(filename = cache_filename, size_limit = size_limit)
                #cache_filenames[cache_filename].restore()
                recache = None
            cache = cache_filenames[cache_filename]

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
        file_header = get_fof_file(group_file,  use_cache = use_cache)
        glen_file = file_header.header.npart[2] #len(ids_in_file)
        if debug>0:
            print( '# get_halo_ids:' ,group_file, 'halo goff',goff, 'halo glen',glen, 'file cumulative goff', goff_file,' file glen',glen_file)
        if True:
            if (goff+glen<=goff_file):
                if debug>0:
                    print( '# get_halo_ids:   (goff+glen<=goff_file) => we read all IDs')
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
                ids_in_file = read_new(group_file,  'PID ', 2, use_cache = use_cache)
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
            if (grnr[isubhalo_in_file]<ihalo):
                continue
            elif (grnr[isubhalo_in_file]>ihalo):
                break
            isubhalo+=1
            subhalo['GRNR'] = ihalo
            subhalo['ihalo'] = ihalo
            subhalo['ifile'] = ifile
            subhalo['isubhalo'] = isubhalo
            if with_ids:
                subhalo['ids'] = halo_ids[subhalo['SOFF']-halo_goff:subhalo['SOFF']+subhalo['SLEN']-halo_goff]
            yield subhalo


def find_progenitors_of_halo(halo, groupbase_format, snap_to=0, max_distance=500., trial=None, max_trials=0, blocks=None, ids_min_frac=0.5, snap_from = None, min_mass=None, use_cache=False, debug=False):
    """
    This function recursively search the progenitor halo in a previous snapshots.
    To do so it loops over all haloes of the new snapshot until it finds one that share a large fraction of particle IDs.
    Once a progenitor is found, the routine search the pro-progenitor in the previous-previous snapshot.
    parameters:
    - `halo` (encoded in a dict,  as returned from yield_halo(...),
    - groupbase_format a routine that, when called for instance like groupbase_format('061') will return the path of the '061' catalog.
    - final snapshot where to search for progenitor (in doubt set to 0)
    - use_cache, set True or to a string to speed up the process
    """
    halo_pos = halo['GPOS']
    halo_ids = halo['ids']

    if trial is None:
        trial = max_trials
    
    # potrei volerlo passare durante i tentativi falliti
    if snap_from is None: 
        snap_from = halo['snap_from']
 
    # we reached the final snapshot
    if (snap_from-1)<=snap_to:
        return

    snap_now = snap_from-1
    groupbase2 = groupbase_format(snap_now)

    if debug:
        print('# I search in snapnum:', snap_now, ' haloes with IDs from prev. halo')

    for halo2 in  yield_haloes(groupbase2, with_ids=True, blocks=blocks, use_cache= use_cache):
        
        # non ha senso andare piu giu
        if halo2['GLEN']<1e4: 
            break
        
        ids_len = len(np.intersect1d(halo_ids, halo2['ids'], assume_unique=True))
        ids_frac = float(ids_len)/float(len(halo_ids))
        if ids_frac>ids_min_frac:
            halo2['snap_from'] = snap_now
            halo2['ids_frac'] = ids_frac

            if debug:
                print('# found match: ihalo: ', halo2['ihalo'],', snapnum:', snap_now, ', perc. of IDs in common: %.2f'%ids_frac)
            yield halo2
            
        
            # abbiamo finito con la ricerca della MAH
            if(halo['MCRI']<min_mass):
                if debug:
                    print('# we found the least massive progenitor')
                return
            
            yield from find_progenitors_of_halo(halo2, groupbase_format,  snap_to,
                                                max_distance=max_distance, max_trials = max_trials,
                                                min_mass = min_mass,
                                                blocks=blocks,
                                                ids_min_frac=ids_min_frac, debug=debug, use_cache = use_cache)
    if debug:
        print('# we didnt find any match in snapnum', snap_now, ' has %.2f'%ids_frac, '% of common IDs, we will try on ',trial,'on previous catalogs')
    # we didnt find any halo, let's try our luck in the next timeslice
    if trial>0:
        if debug:
            print('# we will search for the same IDs in snapnum:', snap_now-1)
        yield from find_progenitors_of_halo(halo, groupbase_format,  snap_to,
                                        min_mass = min_mass,
                                            max_distance=max_distance*1.2, trial=trial-1, max_trials = max_trials, 
                                            blocks=blocks, snap_from = snap_now,      ids_min_frac=ids_min_frac, debug=debug, use_cache = use_cache)
    else: # we give up if we dont find any halo in `_trial_default` snapshots back
        if debug:
            print('# no more progenitors found')
        return

@memoize
def match_from_halo(halo_pos, halo_ids, groupbase2,   blocks = None,  max_distance=500., ids_min_frac=0.3, min_glen=1e4, use_cache = False):
    """
    given a position and a list of IDs, search a halo in a given groupbase (groupbase2). This routine is used by yield_matches_from_haloes
    to search for the same halo in two different simulated snapshot that were run with the same ICs and possibily with different physics
    """
    if blocks is None:
        blocks = ('GPOS','GLEN')

    for halo2 in  yield_haloes(groupbase2, with_ids=True, blocks=blocks, use_cache = use_cache):

        if halo2['GLEN']<min_glen: # non ha senso andare piu giu
            break

        ids_len = len(np.intersect1d(halo_ids, halo2['ids'], assume_unique = True))
        ids_frac = ids_len/len(halo_ids)
        if ids_frac>ids_min_frac:
            halo2['ids_frac'] = ids_frac
            return halo2




def yield_matches_from_haloes(groupbase1, groupbase2,  blocks = None, filter = None, max_distance=500., ids_min_frac=0.3, min_glen=1e4, use_cache = False):
    """
    Search for the same halo in two different simulated snapshot that were run with the same ICs and possibily with different physics,
    or between to snapshot catalogs of the same simulation.
    The two catalogs are supposed to be groupbase1 and groupbase2.
    This routine, will return a generator that gives all couple of halo1, halo2 (as returned by the same format of yield_haloes(...) ),
    where halo1 belongs to groupbase1 and halo2 belongs to groupbase2.
    """
    if blocks is None:
        blocks = ('GPOS','GLEN')

    for halo1 in  yield_haloes(groupbase1, with_ids=True, blocks=blocks, use_cache = use_cache):
        if halo1['GLEN']<min_glen: #non ha senso andare piu giu
            break
        if filter is None or filter(halo1):
            halo1_pos = halo1['GPOS']
            halo1_ids = halo1['ids']
            halo2 = match_from_halo(halo1_pos, halo1_ids, groupbase2, blocks, max_distance, ids_min_frac, min_glen, use_cache=use_cache)
            yield halo1, halo2



@memoize
def get_all_files(groupbase, use_cache=False):
    """  given a catalog base path `groupbase`, e.g. 'groups_061/sub_061', return a list of all catalog files.
    it could be for instance  ['groups_061/sub_061'] (if it is a single-file catalog) or ['groups_061/sub_061.0', 'groups_061/sub_061.1'] 
    (if it is a catalog stored in two files)"""
    return list(g3.yield_all_files(groupbase))

    

def get_most_massive_halo(groupbase, iend=10, mblock = 'MVIR', blocks=None, use_cache = False):
    """
    get the most massive halo (as returned by `yield_haloes`) in a catalog stored in the path `groupbase` (you can set `groupbase='sub_061'` even if it is stored in multiple files)
    based on the block `mblock` (default value 'MVIR').
    """
    hmax = None
    if blocks is None:
        blocks = (mblock,)
    for _h in yield_haloes(groupbase, 0,iend, blocks=blocks, use_cache=use_cache):
        if hmax is None or _h[mblock]>hmax[mblock]:
            hmax = _h
    return hmax


def get_bcg(groupbase, halo, blocks, use_cache = False):
    """
    given the base path of a catalog  (e.g., you can set `groupbase='sub_061'` even if it is stored in multiple files),
    and a `halo` dict (as returned by `yield_haloes(...)`), it returns its most massive subhalo (as returned by `yield_subhaloes`)
    """ 
    for sh in yield_subhaloes(groupbase,  ihalo=halo['ihalo'],
                                blocks = blocks,
                                halo_goff = halo['GOFF'],
                                use_cache=use_cache):
            return sh
    
                

"""

###############################################################################
#            FROM HERE BELOW A LIST OF UNDOCUMENTED FUNCTIONS                 # 
###############################################################################

"""

@memoize
def read_particles_in_box(snapbase, gpos, radius, blocks, ptypes, use_super_indexes = None, use_cache=False):
    return g3.read_particles_in_box(snapbase, gpos, radius, blocks, ptypes, use_super_indexes = use_super_indexes)

def get_subhalo_particles(groupbase,  subhalo, halo, p = None, cache_p = None, use_cache = None):
    soff =  subhalo['SOFF']
    slen =  subhalo['SLEN']

    grnr = subhalo['ihalo']
    sub_id = subhalo['isubhalo'], halo['ihalo']
    subhalo['_'] = subhalo['isubhalo'],halo['ihalo']

    if cache_p is not None:
            
        _cache = cache_p
        
    else:
        _cache = {}                
    if 'nearby_haloes' not in _cache:
        _cache['nearby_haloes']  = {}
    if 'subs' not in _cache:
        _cache['subs']  = {}
   
    if p is None:
        raise Excepiton("to run get_subhalo_particles you must first run get_halo_particles and either use the same use_cache or pass its result in the p parameter")
        
    if  grnr not in _cache['nearby_haloes']:

        _halo = next( yield_haloes(groupbase, grnr, with_ids=False, use_cache = use_cache))
        if _halo['ihalo']!=grnr:
            raise Exception('ihalo != grnr, %d,%d'%(halo['ihalo'], grnr))

        _cache['nearby_haloes'][grnr] ={
            'ids': get_halo_ids(groupbase, halo['GOFF'], halo['GLEN'], use_cache = use_cache)[0],
            'goff':  halo['GOFF'],
            'glen': halo['GLEN']
        }
        if cache_p is not None:
            _cache.store()


    data_p = _cache['subs']
    if sub_id not in data_p:
        nearby_halo_ids = _cache['nearby_haloes'][grnr]['ids']
        nearby_halo_goff = _cache['nearby_haloes'][grnr]['goff']
        nearby_halo_glen = _cache['nearby_haloes'][grnr]['glen']
        particelle_all =  p
        


        subhalo_ids =  nearby_halo_ids[int(soff - nearby_halo_goff):int(soff + slen - nearby_halo_goff)]

        subhalo_particle_masks = np.in1d(particelle_all['ID  '],subhalo_ids)


        data_p[sub_id] = {"p": {k:  particelle_all[k][subhalo_particle_masks] for k in ['POS ','MASS','PTYPE']}}

        if cache_p is not None:
            _cache.store()
                                    

    return data_p[sub_id]['p']


def set_p(d, k, f, *l, v=0, on_edited = None, debug = False, **kw):
    d.setdefault('_edited', False)
    d.setdefault('_versions', {})
    _v = d['_versions'].setdefault(k, None)
    if v==_v:
        return d[k] 
    if callable(f):
        try:
                d[k] = f(*l, **kw)
        except:
            print(v)
            print(l)
            print(dict(kw))
            raise
    else:
        d[k] = f
    d['_versions'][k] = v
    d['_edited'] = True
    if on_edited is not None:
        if debug:
            print('# store in set_p',k)
        on_edited()
    return d[k]

def on_edited(d, f, debug=False, layer = 0):
    run = False
    for k,v in d.items():
        if k=='_edited' and v:
            if debug:
                print("# set run in ",k, '_edited:',v)
            d['_edited'] = False            
            run = True
        elif isinstance(v, dict):

            _run = on_edited(v, f,  layer = layer+1, debug = debug)

            run = run or _run
            if debug and _run:
                print("# edited in ",k, run)
    if layer == 0 and run:
        if debug:
            print('# exec f..')
        f()
        if debug:
            print('# execd f.')
    return run



def get_subhalo(grouppath,   ihalo, isubhalo, blocks = None, use_cache = None):
    halo = None
    halog = yield_haloes(grouppath, ihalo, ihalo_end=None, with_ids = True, use_cache=use_cache)
    for halo in halog:
        break
    if blocks is None:
        blocks = ('RHMS','SPOS','GRNR','VMAX','MSUB','SVEL')    
    if halo is None: #there are no haloes
        return
    
    subhaloesg = yield_subhaloes(grouppath, ihalo=halo['ihalo'], with_ids = True, halo_ids = halo['ids'],  halo_goff = halo['GOFF'],  
                                          blocks=blocks, use_cache=use_cache) 
    for _isubhalo, subhalo in enumerate(subhaloesg):
        if _isubhalo == isubhalo:
            return subhalo
    raise Exception('%s has no (%d, %d)'%(grouppath,   ihalo, isubhalo   ))


def find_subhalo_in_snap(grouppath,   ids, threshold=0.5, use_cache = None):
    candidates = {}
    print('#        looping ',grouppath)
    for halo  in  yield_haloes(grouppath, 0,  with_ids = True,  ihalo_end=None,
                                      blocks=('GLEN', 'MVIR', 'RVIR', 'GPOS'), use_cache=use_cache):
            for subhalo in  yield_subhaloes(grouppath, ihalo=halo['ihalo'],   with_ids = True, halo_ids = halo['ids'],  halo_goff = halo['GOFF'],  
                                                    blocks=('RHMS','SPOS','GRNR','VMAX','MSUB','SVEL'), use_cache=use_cache):
                common = np.sum(np.in1d(ids, subhalo['ids']))
                f = float(common)/float(len(ids))
                if f>threshold:
                    print('#            found with ',f)
                    return subhalo
                elif f>0:
                    print('#            store one with ',f)
                    candidates[f] = subhalo
    if len(candidates)>0:
        min_key = max(candidates.keys())
        print('#            return ',min_key)
        return candidates[min_key]
    else:
        return None

def get_subhalo_history(groupformat, snap_start, ihalo, isubhalo, history, use_cache = None):
    subhalo = None
    print('# history of ',groupformat(snap_start),(ihalo,isubhalo))
    for snap_i in  reversed(range(1, snap_start+1)):
        if snap_i in history:
            subhalo = history[snap_i]
            continue
        
        grouppath = groupformat(snap_i)
        print('#    snap',snap_i, grouppath)
        if not os.path.isfile(g3.get_one_file(grouppath)):
            print('#    not exists.. continue')
            continue
        #print(subhalo)
        if subhalo is None:
            subhalo = get_subhalo(grouppath, ihalo, isubhalo, use_cache = use_cache)
        else:
            subhalo = find_subhalo_in_snap(grouppath, subhalo['ids'], use_cache = use_cache)
        if subhalo is None:
            return
        yield snap_i, subhalo

