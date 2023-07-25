
This repository hosts a  collection of tools to read and post-process  large `Gadget2` and `Gadget3` files (including key files).
The core routines (`g3read.py`) are a composition compbining some `pynbody` and a porting of some legacy _Klaus' IDL_ routines.
to send batch jobs to the [c2pap web portal](http://c2papcosmosim.uc.lrz.de/)and to convert gadget files to HDF5

**For questions**: Antonio Ragagnin <antonio.ragagnin@unibo.it> (https://aragagnin.github.io)


**Table of Contents:**

- [Install](#install)
- [read Gadget and key files with g3read.py](#read-gadget-and-key-files-with-g3readpy)
  - [Read a single Gadget file](#read-a-single-gadget-file)
  - [Access the header](#access-the-header)
  - [Reading FOF or Subfind files](#reading-fof-or-subfind-files)
  - [Reading from a large run with super indexes](#reading-from-a-large-run-with-super-indexes)
  - [Writing a new snapshot file from scratch](#writing-a-new-snapshot-file-from-scratch)
- [Looping through haloes, their subhaloes, and IDs with g3matcha.py](#looping-through-haloes-their-subhaloes-and-ids-with-g3matchapy)
  - [Looping through haloes and sub haloes](#looping-through-haloes-and-sub-haloes)
  - [Caching of data to speedup SubFind or FoF reading](#caching-of-data-to-speedup-subfind-or-fof-reading)
  - [ Finding main progenitors](#finding-main-progenitors)
- [batch jobs for http://c2papcosmosim.uc.lrz.de/ with c2pap_batch.py](#batch-jobs-for-httpc2papcosmosimuclrzde-with-c2pap_batchpy)


# Install

The spirit of this collection of files is that you do not need to download all those files: you can read the documentation and just download the single files that you need for your tasks.

However on `python3` you can also install the package with the following command:

```bash
python -mpip install  git+https://github.com/aragagnin/g3read
```

and you will be able to `import g3read`, `import g3matcha` or `import g3read_units` in your scripts without any further download.

# Read Gadget and key files with g3read.py

To read  napshots and FoF/SubFind outputs all you need is `g3read.py`. This library contains the GadgetFile class from [pynbody](https://github.com/pynbody/pynbody). The library `g3read` will use [numba](http://numba.pydata.org) if available.

## Read a single Gadget file

The easiset way to read from a Gadget file is to use the function `read_new` (a clone of Klaus Dolag IDL routine). 

```python
read_new(filename, blocks, ptypes, center=None, is_snap=True, join_ptypes=False)
```

- `filename`: path of the gadget file
- `blocks`: block(s) to read. Can be a string or a list of one or more strings (e.g. `"MASS"`, `["MASS"]`, or `["MASS", "POS "]`). Blocks must be 4characters long. E.g. use `"ID   "` to read IDs.
- `ptypes`: can be an integer  or a list of integers representing the particle type (e.g. `1`, `[1]`, or `[0,1,2,3,4,5]`). Using `-1` equals to asking for all blocks.
`center`: if set turns on the periodic (PBC) assumpions. 
- `is_snap`: default is `True`. Set to `False` in order to read SUBFIND data
- `join_ptypes`: default is `True`. If `False` return a dictionary with data separated by ptype.
- return type depends on the input data: if `blocks` is a list, then the result is a dictionary of array data for each block (see examples below) 

Examples:

```python
import g3read

# read one block per time from all ptype:
mass =  g3read.read_new("./test/snap_132", "MASS", -1) #the -1 means to select all particles
pos  = g3read.read_new("./test/snap_132", "POS ", -1) 
x = pos[:,0]
y = pos[:,1]

# read one block per time from one ptype:
mass =  g3read.read_new("./test/snap_132", "MASS", 4) #the 4 means to select star particles
pos  = g3read.read_new("./test/snap_132", "POS ", 4) 

# read multiple blocks from one ptype
data =  g3read.read_new("./test/snap_132", ["POS ", "MASS"], 0) #the 0 select only gas particles
pos = data["POS "]
mass  = data["MASS"]

# read multiple blocks from multiple ptypes, all stacked in a single array
data =  g3read.read_new("./test/snap_132", ["POS ", "MASS","TEMP"], [0,4]) #the 0,4 select  gas and star particles
stars_and_gas_pos = data["POS "]
stars_and_gas_mass = data["POS "]
#And use the artificial block `PTYPE` to filter by particle type:
gas_mask = data["PTYPE"]==0
gas_temp = data["TEMP"][gas_mask] 

# read multiple blocks from multiple ptypes, each in a separate array (note join_ptyes=False)
data =  g3read.read_new("./test/snap_132", ["POS ", "MASS"], [0,4], join_ptypes=False) #the 0,4 select  gas and star particles
stars_pos = data[4]["POS "]
stars_mass  = data[4]["MASS"]
gas_pos = data[0]["POS "]
gas_mass  = data[0]["MASS"]

```

There is **Object Oriented version** in case you need multiple call `read_new` multiple times from a file, you can instantiate a GadgetFile  separately (note `is_snap=False` when reading a catalog):

```python
f = g3read.GadgetFile("./test/snap_132")
mass =  f.read_new("MASS", 4)
pos  = f.read_new("POS ", 4) 
```
## Access the header

You can access the header field of the object `GadgetFile`. As in the following example

```python
f = g3read.GadgetFile("./test/snap_132")
numpart_thisfile = f.header.npart
numpart_allfiles = f.header.npartTotal
redshift =  f.header.redshift
scale_factor =  f.header.time
box_size =  f.header.BoxSize
h0 = f.header.HubbleParam
all_other_properties_name = f.header.__dict__.keys()
```

## Reading FOF or Subfind files

To read from the catalog you need to use `read_new` with the flag `is_snap=False`.
When reading this documentation please pay attention that there is an ambiguity when mentioning Subfind:
it is both the algorithm to find only subhaloes within a halo (as opposed to FoF that only searches for haloes), 
however the version implemented in gadget does write catalog files that contain also information about the FoF groups.
There fore here below we will read both FoF haloes and Subfind subhaloes by reading Subfind outputs.

```python
import g3read
snapbase = '/HydroSims/Magneticum/Box2/hr_bao/snapdir_136/snap_136'
groupbase = '/HydroSims/Magneticum/Box2/hr_bao/groups_136/sub_136'

# note that here set ptype as 0  (third argument)  because Subfind output of the FoF catalog
# put it in the zeroth particle block. You would have to do the same even if you read the output
# using Idl routines
fof =  g3read.read_new(groupbase+'.0', ['GPOS','RCRI','MCRI'], 0, is_snap=False) 

# here we are reading the subhaloes from Subfind catalog
subfind =  g3read.read_new(groupbase+'.0', ['SPOS','SVEL','GRNR'], 1, is_snap=False) 

print('positions: ', fof['GPOS'])
print('r200cs: ', fof['RCRI'])
print('m200cs ', fof['MCRI'])

```

Check `test_g3read.py` for a sample that converts subfind haloes to ASCII table.

Sometimes your simulation produce a small number of catalog files, g3read can concatenate the results for you if you use the parameter `multiple_files`:

```python
fof =  g3read.read_new(groupbase, ['GPOS','RCRI','MCRI'], 0, is_snap=False, multiple_files=True ) 
subfind =  g3read.read_new(groupbase, ['SPOS','SVEL','GRNR'], 1, is_snap=False, multiple_files=True )) 

```




## Reading from a large run with super indexes

The signature of `g3read.read_particles_in_box` is almost the same of `read_new`. As opposed to `read_new`, `g3read.read_particles_in_box` additionally needs a minimum radius.

```python
read_particles_in_box(snap_basepath, center, distance, blocks, ptypes, join_ptypes=True, is_snap=True)
```

- `snap_basepath`: is the base path of the snapshot, (e.g. `/HydroSims/Magneticum/Box2/hr_bao/snapdir_136/snap_136` )
- `center`: position of the sphere of radius `distance` that will contains our particles
- `radius`: searching radius  of particles. Note the routine will return a superset of particiles within `radius`.
- `blocks`: as `read_new`
- `ptypes`: as `read_new`
- `is_snap`:  as `read_new`
- `join_ptypes`:  as `read_new`
- return  as `read_new`

In this example I first read the position and radius of a FoF object (from the fof files) and then I extract its properties with `read_particles_in_box`.

```python
import g3read
snapbase = '/HydroSims/Magneticum/Box2/hr_bao/snapdir_136/snap_136'
groupbase = '/HydroSims/Magneticum/Box2/hr_bao/groups_136/sub_136'

#read FOF data
fof =  g3read.GadgetFile(groupbase+'.0', is_snap=False) #if you read a FoF/Subfind file, add is_snap = False 
halo_positions = fof.read("GPOS",0) #block zero has FoF data, block 1 has SubFind data
halo_radii = fof.read("RVIR",0)

#extract position of first halo
first_halo_position = halo_positions[0]
first_halo_radius = halo_radii[0]

#use read_particles_in_box
f = g3read.read_particles_in_box(snapbase, first_halo_position, first_halo_radius, ["POS ","MASS"], -1)
x = f["POS "][:,0]
y = f["POS "][:,1]
mass = f["MASS"]

group_gas_data = g3.read_particles_in_box(snapbase, first_halo_position, first_halo_radius,   ['MASS', 'TEMP', 'POS '], 0)

#note: in the spirit of Klaus read_particles_in_box, the above routine returns a superset of particles within `r200c`.
#we now filter data outside r200c, we use  g3.to_spherical that returns an array of [rho, theta, phi] around `center`.

group_gas_distance_from_center = g3.to_spherical(group_gas_data['POS '], center).T[0]
group_gas_mask = group_gas_distance_from_center < r200c
group_gas_masswtemp =  group_gas_data['TEMP'][group_gas_mask] * group_gas_data['MASS'][group_gas_mask]
group_gas_avgtemp =  np.mean(group_gas_masswtemp)/np.sum( group_gas_data['MASS'][group_gas_mask])
```

Note you could also call `fof =  read_new(groupbase, "GPOS", 0, multiple_files=True, is_snap=False)` to join all FoF catalog files (`sub_136.0`, `sub_136.1`, etc..). You should do it only if you know it can fit your memory.

## Writing a new snapshot file from scratch

The class `g3read.GadgetFile` has a function `write_block` that will overwrite a block with a new provided array.

In this example I recompute the gravitational potentail between particles and store it back in a pre-existing `"POT "` block.
The pacakge `pp.py` contains a routine to compute the gravitational potential between particles of the snapshot.


```python
import g3read, pp 
my_filename = "./test/snap_132"
my_filename_output "./test/new_snap_132"

f = g3read.GadgetFile(my_filename)
positions = f.read_new("POS ",-1) #-1 means all particles
masses = f.read_new("MASS",-1)
potential = pp.gravitational_potential(masses, positions, center).potential

f.write_block("POT ", -1, potential, filename=my_filename_output)
```


Here an example on how to create a Gadget snapshot and its header from scratch

```python

import g3read as g3, numpy as np

### START INPUT DATA ###

filename = 'mysnap'

npart = np.array([10,20,0,0,0,0])
mass_table  =  [0.]*6
mass_table[1] = 0.1 #ptype1 has fixed mass == 0.1
redshift = 1.0
time = 1./(redshift+1.)
BoxSize=100.
Omega0 = 0.27
OmegaLambda = 1. - Omega0
HubbleParam = .704 #or should I put 70.4? I do not remember
num_files = 1

### END INPUT DATA ###

with open(filename, 'a') as f: #create file if doesn't exists
    pass

# generate header
header = g3.GadgetHeader(npart, mass_table, time, redshift, BoxSize, Omega0, OmegaLambda, HubbleParam, num_files=num_files)
# write header to file
f = g3.GadgetWriteFile(filename, npart, {}, header) #write header file
f.write_header(f.header)

#allocate blocks data
f.add_file_block('POS ', 30*4*3  , partlen=4*3) #add a block of 30*4*3 bytes each of 4*3 bytes
f.add_file_block('MASS', 10*4, partlen=4) #add a block of 10*4*3 bytes each of 4 bytes

#write blocks data to disk
f.write_block( 'POS ', -1, np.array([[1,2,3]]*30, dtype=np.float32)) #write 30 positions
f.write_block( 'MASS', -1, np.array([1,2,3,4,5,6,7,7,7,7], dtype=np.float32)) #write 10 masses
```
# Looping through haloes, their subhaloes, and IDs with g3matcha.py

`g3matcha.py` (which depends on `g3read`) provides high level functionality to perform for loop over haloes and their subhaloes.
Check file `test_g3read_ids.py` for a complete example.

## Looping through haloes and sub haloes

To loop over haloes use `yield_haloes`, with the following parameters:

- `group_base`: the base path of the SubFind group file, e.g. `groups_000/sub_000`
- `ihalo_start`: first halo index to be returned (default 0)
- `ihalo_end`: maximum halo index to be read (default None, to read all haloes)
- `min_mcri`: skip all haloes with a mass lower than `min_mcri`
- `blocks`: FoF blocks to be read, default: `["GPOS", "RCRI", "MCRI", "GOFF", "GLEN"]`
- `with_ids`: returns also the halo IDs in the block `"ids"` (its lowercase and of three character to stress that it is artificial), default is `False`.

and as result yields a list of haloes with the chosen properties plus the index of halo `ihalo`, the SubFind file index that contains it `ihalo_file`, its position in said file `ihalo_in_file` (start from 0 at each SubFind file) and `boxsize` from the file header (see previous section).

To loop over haloes use `yield_subhaloes`, with the following parameters:

- `group_base`: as `yield_haloes`
- `ihalo`: index of parent `halo`
- `ifile_start`: from which SubFind partial file it must start reading (typycally you can set it as `parent_halo['ihalo_file']` to speed computation), default is 0
- `whith_ids`: if true, each returned subhalo will contian also their IDs. You must provide the parent halo IDs with the keyword `halo_ids`. IDs will be provided with tbe aritificial block `ids`
- `halo_ids`: halo ids of the parent halo

It returns a list of subhaloes of the given parent haloes. 

For instance, let's expand the previous example with the reading of subhaloes:

```python
import g3read as g3, g3matcha as matcha, numpy as np
groupbase = '/HydroSims/Magneticum/Box1a/mr_bao/groups_144/sub_144'
# read first 10 fof haloes
for halo  in  matcha.yield_haloes(groupbase, 0, ihalo_end, blocks=('GLEN', 'MCRI', 'RCRI', 'GPOS')):
    print('halo number:', halo['ihalo'], halo)
    #now loop over subhaloes of the parent halo
    for subhalo in  matcha.yield_subhaloes(groupbase, ihalo=halo['ihalo']):
        print('    sub halo information:', subhalo)
```    

Here below an example that read both FoF and subhalo IDs. Note that we must provide the infromation ` halo_goff = halo['GOFF']` to `yield_subhaloes` in order to grab sub halo IDs.

```python
import g3read as g3, g3matcha as matcha, numpy as np
groupbase = '/HydroSims/Magneticum/Box1a/mr_bao/groups_144/sub_144'
# read first 10 fof haloes
for halo  in  matcha.yield_haloes(groupbase, 0, with_ids=True, ihalo_end=10, blocks=('GLEN', 'MCRI', 'RCRI', 'GPOS')):
    print('halo number:', halo['ihalo'])
    print('halo IDs: ', halo['ids'])
    #now loop over subhaloes of the parent halo
    for subhalo in  matcha.yield_subhaloes(groupbase,  with_ids=True, halo_ids = halo['ids'],  halo_goff = halo['GOFF'], ihalo=halo['ihalo']):
        print('    sub halo IDs:', subhalo['ids'])
```

## Caching of data to speedup SubFind or FoF reading

`g3matcha` routines can be cached  to a dict in order to make reading  by adding `use_cache=True` to function calls.
You can also cache results to file in order to recycle reads when running the same script multiple time (we all know you'll run your script many many times) by providing  provide a `g3matcha.cache_filename`.

Here how it will look like in the previous example

```python
import g3read as g3, g3matcha as matcha, numpy as np

#set to None if you do not want to cache data to file (cache will stay in memory then avaiable only for this run)
matcha.cache_filename = 'cache'

groupbase = '/HydroSims/Magneticum/Box1a/mr_bao/groups_144/sub_144'
# read first 10 fof haloes
for halo  in  matcha.yield_haloes(groupbase, with_ids=True, ihalo_end=10, blocks=('GLEN', 'MCRI', 'RCRI', 'GPOS'), use_cache= True):
    print('halo number:', halo['ihalo'])
    print('halo IDs: ', halo['ids'])
    #now loop over subhaloes of the parent halo
    for subhalo in  matcha.yield_subhaloes(groupbase,  with_ids=True, halo_ids = halo['ids'],  halo_goff = halo['GOFF'], ihalo=halo['ihalo'], use_cache = True):
        print('    sub halo IDs:', subhalo['ids'])

# the script will be much faster now!
```

## Finding main progenitors

The example below will find all progenitors of the main halo at snapshot 091.

```python

#we let g3matcha.find_progenitors_of_halo know how to generate FoF catalog path 
groupbase_format = lambda snap_num: './groups_%03d/sub_%03d'%(snap_num, snap_num)

snap_from = 91 # snapshot number where to find the target halo (progenitors will start from snap 62)
halo_nr = 0 # index of target halo in snap 063
snap_to=20 #lower snapshot number where to find for progenitors

blocks = ('GPOS','RVIR','GLEN','MVIR') #blocks to be displayed per halo
    
# if one progenitor is not found in the previous snapshots, then code will keep search on `_trial_default` snapshots before it. 
_trial_default = 3

min_mass = 1e3. #minimum mass of the oldest progenitor

groupbase = groupbase_format(snap_from)
halo = next(g3m.yield_haloes(groupbase,  ihalo_start = halo_nr,  blocks=blocks, with_ids = True))
halo['snap_from'] = snap_from

print('# Starting loop on progenitors: ')

for progenitor in g3matcha.find_progenitors_of_halo(
        halo,
        groupbase_format,
        snap_to,
        use_cache = False, #set to true to speedup
        max_trials = _trial_default,
        blocks=blocks,
        min_mass = min_mass):
        
    print("- ihalo: ", progenitor['ihalo'])
    print("  snap: ", progenitor['snap_from'])
    print("  ids_frac: ", progenitor['ids_frac'])
```

#  batch jobs for http://c2papcosmosim.uc.lrz.de/ with c2pap_batch.py

Given a list of clusters previously extracted from the c2pap web portal (the output file name is `dataset.csv`), the script `c2pap_batch.py` automatize the process of sending the same jobs parameter to all those haloes.


The script takes the parameter `-s <Service name>` where the service name can be `SMAC`,`SimCut`,`PHOX`,`query`,
the path of the dataset file name `-f <daraset.csv>`, the job parameters must be set with `-p` and must
consist of the value displayes in the web portal.

Below a list of all parameters. **Note: they are case sensitive.**

```
                        Form parameters.
                            SimCut:
                            IMG_Z_SIZE
                            r500factor

                            SMAC:
                            content
                            IMG_SIZE
                            IMG_Z_SIZE
                            PROJECT
                            r500factor

                            PHOX:
                            mode
                            instrument (-1 for generic)
                            instrument_a (only if generic)
                            instrument_fov (only if generic)
                            t_obs_input
                            img_z_size
                            simulate
```

For instance the following will run a SMAC job over all objects in the dataset.csv:

```bash
python c2pap_batch.py -f dataset.csv -s SMAC -p content="bolometric x-ray luminosity" IMG_SIZE=512 IMG_Z_SIZE=5000 PROJECT='along z, xy plane' r500factor=2.0 -u <YOUR USERNAME> 
````

To avoid running duplicate jobs, use the flags ` --existing-jobs --cache-jobs cachefile.pickle` to make the script check for existing jobs with identical parameters.

