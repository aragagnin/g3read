"""
Gadget high level reader functions with units.

"""

from pint import Context
from pint import UnitRegistry
import g3read as g3



class Units(object):
    def __init__(self):
        self.conversion_blocks = {
            'MASS':'gmass',
            'POS ': 'glength',
            'RHO ': 'gmass/glength^3'
        }
        self.ureg = UnitRegistry()
        u = self.ureg
        u.define('Msun = 1.99885e30kg')
        u.define('gmass = 1e10 Msun/hubble')
        u.define('cmass = Msun/hubble')
        u.define('clength = kpc/hubble*scalefactor')
        u.define('glength = clength')
        u.define('cvelocity = scalefactor*km/s')
        u.define('gvelocity_a = (scalefactor**0.5)km/s')
        u.define('gvelocity_noa = km/s')

        #u.add_context(c)
    def set_context(self, hubble, scalefactor, debug=False):
        u = self.ureg
        #c = self.context
        #u.enable_contexts(c,hubble=hubble, scalefactor=scalefactor)
        if(debug):
            print('# hubble = ',hubble)
            print('# scale factor = ',scalefactor);
            print('# redshift = ',1./scalefactor-1.);
        u.define('hubble = '+str(hubble))
        u.define('scalefactor = '+str(scalefactor))
    def get_u(self):
        return self.ureg

def gen_units(f,units=None, debug=False):
    if units is None:
        units = Units();
        units.set_context(f.header.HubbleParam, f.header.time, debug=True);
    return units


def gen_factor(units,f, blocks):
    units = gen_units(f,units)
    ureg = units.get_u();

    factor = {}
    for block in g3.iterate(blocks):
        factor[block] = ureg.parse_expression(units.conversion_blocks[block])
    return factor

def get_units(snap_path):
    one_file_name = g3.get_one_file(snap_path)
    one_file_descriptor = g3.GadgetFile(one_file_name)
    units = gen_units(one_file_descriptor, debug=True)
    return units


def read_new(filename, blocks, ptypes, join_ptypes=True, only_joined_ptypes=True, periodic=True, center=None, is_snap=False, units=None):
    f = g3read.GadgetFile(filename)
    factor = gen_factor(units, f, blocks)
    res = read_new(filename, blocks, ptypes, join_ptypes=join_ptypes, only_joined_ptypes=only_joined_ptypes, periodic=periodic, center=center, is_snap=is_snap, factor=factor)
            

def yield_particles_blocks_in_box(snap_file_name,center,d, blocks, ptypes,  periodic=True, units=None, debug=False):
    #we get the first filename to get aGadgetFile object and get hubble/scalefactor
    if debug: print('# wait for keyfiles of ',snap_file_name,center)
    filename = None
    keylist = None
    part_keylist = [x for x in g3.yield_particles_file_in_box(snap_file_name,center,d,debug=debug)]
    filename = part_keylist[0][0] #first 0 is for the first, pair file-keylit, second zero is to grab the file
    if debug: print('# found ',filename)
    f = g3.GadgetFile(filename)
    factor = gen_factor(units, f, blocks)
    if debug: print('# factor',factor)
    yield from g3.yield_particles_blocks_in_box(snap_file_name,center,d, blocks, ptypes,  periodic=periodic, factor=factor, debug=debug)# part_keylist = part_keylist)
     

def read_particles_in_box(snap_file_name,center,d, blocks, ptypes, join_ptypes=True, only_joined_ptypes=True, periodic=True):
    for filename in g3.yield_particles_file_in_box(snap_file_name,center,d) :
        break
    f = g3read.GadgetFile(filename)
    factor = gen_factor(units, f, blocks)
    return  g3.yield_particles_blocks_in_box(snap_file_name,center,d, blocks, ptypes,  periodic=periodic, factor=factor)
