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
                cache_filenames[cache_filename].restore()
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
            subhalo['GRNR'] = ihalo
            subhalo['ihalo'] = ihalo
            subhalo['ifile'] = ifile
            subhalo['isubhalo'] = isubhalo
            if with_ids:
                subhalo['ids'] = halo_ids[subhalo['SOFF']-halo_goff:subhalo['SOFF']+subhalo['SLEN']-halo_goff]
            yield subhalo


def find_progenitors_of_halo(halo, groupbase_format,snap_to, max_distance=500., trial=None, max_trials=0, blocks=None, ids_min_frac=0.5, snap_from = None, min_mass=None, use_cache=False, debug=False):

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
