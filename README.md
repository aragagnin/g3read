# g3read

These tools give you the possibility to read and do post processing of large Gadget2 and Gadget3 files (including key files),  to send batch jobs to the <a href="http://c2papcosmosim.uc.lrz.de/" rel="nofollow">c2pap web portal</a> and to convert gadget files to HDF5.

You do not need to download all those files: this is a collection of libraries, so read the documentation and just download what you need for your task.


Table of Contents
=================

   * [g3read](#g3read)
      * [Read a single Gadget file](#read-a-single-gadget-file)
      * [Writing  back to a (new) file](#writing--back-to-a-new-file)
      * [Reading from a large run (with super indexes)](#reading-from-a-large-run-with-super-indexes)
      * [High Performance](#high-performance)
      * [Working with units of measurements](#working-with-units-of-measurements)
      * [Maps of large simulations](#maps-of-large-simulations)
   * [Submit a batch of jobs to the c2pap web portal (<a href="http://c2papcosmosim.uc.lrz.de/" rel="nofollow">http://c2papcosmosim.uc.lrz.de/</a>)](#submit-a-batch-of-jobs-to-the-c2pap-web-portal-httpc2papcosmosimuclrzde)
   * [Convert Gadget2/3 files to HDF5](#convert-gadget23-files-to-hdf5)
   * [Read Gadget files with units](#read-gadget-files-with-units)
   * [Create SMAC-like maps from large samples](#create-smac-like-maps-from-large-samples)


## Read a single Gadget file

To read  napshots and FoF/SubFind outputs all you need is `g3read.py`. This library contains the GadgetFile class from [pynbody](https://github.com/pynbody/pynbody).

The easiset way to read from a Gadget file is to use the function `read_new` (a clone of Klaus Dolag IDL routine). 

```python
read_new(filename, blocks, ptypes, center=None, is_snap=False)
```

- `filename`: path of the gadget file
- `blocks`: block(s) to read. Can be a string or a list of one or more strings (e.g. `"MASS"`, `["MASS"]`, or `["MASS", "POS "]`)
- `ptypes`: can be an integer  or a list of integers representing the particle type (e.g. `1`, `[1]`, or `[0,1,2,3,4,5]`). Using `-1` euals to asking for all blocks.
`center`: if set turns on the periodic (PBC) assumpions. 
- return type depends on the input data: if `blocks` is a list, then the result is a dictionary of array data for each block (see examples below) 

Example:

```python
import g3read
mass =  g3read.read_new("./test/snap_132", "MASS", -1) #the -1 means to select all particles
pos  = g3read.read_new("./test/snap_132", "POS ", -1) 
x = pos[:,0]
y = pos[:,1]
```

To select only gas particles:

```python
pos_gas =  g3read.read_new("./test/snap_132", "POS ", 0) #the 0 select only gas particles
````


Note the difference if `"POS "` is within a list:

```python
data =  g3read.read_new("./test/snap_132", ["POS "], 0)  
pos = data["POS "]
````

To select both position and mass of gas and dark matter:

```python
data =  g3read.read_new("./test/snap_132", ["POS ", "MASS"], [0,1]) #the 0 select only gas particles
pos = data["POS "]
mass  = data["MASS"]
````


In case you need multiple reads and want to save some I/O time,  , you can instantiate a GadgetFile  separately:

```python
`f = g3read.GadgetFile("./test/snap_132")
pos_gas =  f.read_new( "POS ", 0) #the 0 select only gas particles
[...]
mass_gas =  f.read_new("MASS", 0)
```

Use the block `PTYPE` to filter by particle type:

```python
f = g3read.GadgetFile(filename)
data = f.read_new(blocks=["POS ","VEL ","TEMP","MASS"], ptypes=[0,1,2,3,4,5])
gas_temp = data["TEMP"][data["PTYPE"]==0] #zero is for gas,1 dm, 4 stars, 5 BH
```

## Writing  back to a (new) file

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
## Reading from a large run (with super indexes)

The signature of `g3read.read_particles_in_box` is almost the same of `read_new`.
As opposed to `read_new`, `g3read.read_particles_in_box` additionally needs a minimum radius.

In this example I first read the position and radius of a FoF object (from the fof files) and then I extract its properties with `read_particles_in_box`.

```python
import g3read
snapbase = '/HydroSims/Magneticum/Box2/hr_bao/snapdir_136/snap_136'
groupbase = '/HydroSims/Magneticum/Box2/hr_bao/groups_136/sub_136'
fof =  g.GadgetFile(groupbase+'.0', is_snap=False) #if you read a FoF/Subfind file, add is_snap = False 

halo_positions = fof.read("GPOS",0) #block zero has FoF data, block 1 has SubFind data
halo_radii = fof.read("RVIR",0)

#extract position of first halo
first_halo_position = halo_positions[0]
first_halo_radius = halo_radii[0]

f = g3read.read_particles_in_box(snapbase, first_halo_position, first_halo_radius, ["POS ","MASS"], [0,1,2,3,4,5])
x = f["POS "][:,0]
y = f["POS "][:,1]
mass = f["MASS"]
```

## High Performance 

The library `g3read` will use [numba](http://numba.pydata.org) if available.

## Working with units of measurements

The library `g3units` read GadgetFiles (using `g3read`) and store blocks in [pint](https://pint.readthedocs.io/) datastructures. Gadget length blocks (e.g. `POS `) uses `pint`  units `glength` (defined ad `kpc * scalefactor / hubble`) and masses use `gmass` (degined as  `1e10 Msun/hubble`)

In this example we read code-units data and return it in physical units automatically:

```python
import g3read_units as g3u
snap_base = '/HydroSims/Magneticum/Box2/hr_bao/snapdir_136/snap_136'
units = g3u.get_units(snap_base)
ureg = units.get_u() #pint ureg
center = [500.,500.,500.] * ureg.glength  #we give a center in code units 
distance = 500. * ureg.kpc # distance we give in real kpc

data = g3u.read_particles_in_box(snap_base, center, distance, ["MASS", "POS "], -1):
mask = ((data["POS "][:,0]-center)<distance) & ((data["POS "][:,1]-center)<distance) ((data["POS "][:,2]-center)<distance)   
total_mass = np.sum(data["MASS"])
print(' Total Mass in physical Msun:', total_mass.to('Msun')) 
``` 

## Maps of large simulations

If you use [SMAC](https://wwwmpa.mpa-garching.mpg.de/~kdolag/Smac/) you know that  you can only do 2D maps with a number of particles that fits your RAM memory. `make_maps.py` is slightly compatible with SMAC and is capable of producing maps of objects that do not fit RAM memory.

`make_maps.py` do also smooth gas properties as SMAC. Here below a comparison between the two: 

![comparison between make_maps.py on the left and SMAC on the right](https://i.imgur.com/xmCauqV.png)

I am still not sure why `make_maps.py` do produce a more noisy output. It may be that SMAC uses a different SPH kernel (I use top hat by now).

`make_maps` input parameters do contain units, so there is no more confusion on what kind of data you are passing. Here below the input parameter used to make the image:

```bash

This is the input file used

```bash
IMG_XY_SIZE = 2209.112 glength
IMG_Z_SIZE = 400.0*1.1344 kpc#glength
CENTER_X =   456582.8
CENTER_Y =   220605.1
CENTER_Z =   279066.1
SNAP_PATH = /home/moon/ragagnin/mnt/pr62go/Magneticum/Box2b/hr_bao//snapdir_031/snap_031
PREFIX_OUT = spwze7x7kjx78ar5/povero_
IMG_SIZE = 128
PTYPES = 0,1,2,3,4 #,5
JOB = 2DMAP
MAP_DIVIDE_BY_SURFACE = True
PROPERTY = MASS
RESULT_UNITS = Msun #/cm^2 #Msun
```


# Submit a batch of jobs to the c2pap web portal (http://c2papcosmosim.uc.lrz.de/)

Given a list of clusters previously extracted from the c2pap web portal (the output file name is `dataset.csv`), the script `c2pap_batch.py` automatize the process of sending the same jobs parameter to all those haloes.


The script takes the parameter `-s <Service name>` where the service name can be `SMAC`,`SimCut`,`PHOX`,`query`,
the path of the dataset file name `-f <daraset.csv>`, the job parameters must be set with `-p` and must
consist of the value displayes in the web portal.

Below a list of all parameters.

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

                            query:
                            query
                            page
                            limit
```

For instance the following will run a SMAC job over all objects in the dataset.csv:

```bash
python c2pap_batch.py -f dataset.csv -s SMAC -p content="bolometric x-ray luminosity" IMG_SIZE=512 IMG_Z_SIZE=5000 PROJECT='along z, xy plane' r500factor=2.0 -u <YOUR USERNAME> 
````

To avoid running duplicate jobs, use the flags ` --existing-jobs --cache-jobs cachefile.pickle` to make the script check for existing jobs with identical parameters.


# Convert Gadget2/3 files to HDF5

Use the utility `gadget_to_hdf5.py`.

```bash
python gadget_to_hdf5.py infile outfile
```

In case you need to map names diffrently from the default version, have a look at the source code of `gadget_to_hdf5.py` and edit your own mapping.


# Read Gadget files with units

The library `g3units.py` reads 

# Create SMAC-like maps from large samples
