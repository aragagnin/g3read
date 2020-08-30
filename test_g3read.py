""" 

This file tests g3read library by parsing halo from a simulation catalog
and outputs (1) FoF propertis (mass, radius, center), then it reads halo particles with `read_particles_in_a_box`
and outpus (2) halo weighted temperature.

`snapbase` and `groupbase` are meant to work on dorc1.usm.uni-muenchen.de servers. 
change it by hand to test it in your environment.

Antonio Ragagnin

"""

import g3read as g3, numpy as np


snapbase = '/HydroSims/Magneticum/Box1a/mr_bao/snapdir_144/snap_144'
groupbase = '/HydroSims/Magneticum/Box1a/mr_bao/groups_144/sub_144'


""" BEGIN """


ihalo = -1
for group_filename  in  g3.yield_all_files(groupbase):
    print('# reading ', group_filename)
    
    groups_in_file_data = g3.read_new(group_filename, ['RCRI', 'MCRI', 'GPOS'], 0, is_snap=False)
    groups_in_file = len(groups_in_file_data['MCRI'])
    for i in range(groups_in_file):
        ihalo+=1
        center = groups_in_file_data['GPOS'][i]
        r200c = groups_in_file_data['RCRI'][i]
        M200c = groups_in_file_data['MCRI'][i]
        print(' - ihalo: ',ihalo)
        print('   center: ', center, '[code units]' )
        print('   R200c: ', r200c, '[code units]' )
        print('   M200c: ', M200c, '[code units]' )
        
        group_gas_data = g3.read_particles_in_box(snapbase, center, r200c,    ['MASS', 'TEMP', 'POS '], 0)
        
        #note: in the spirit of Klaus read_particles_in_box, the above routine returns a superset of particles within `r200c`.
        #we now filter data outside r200c, we use  g3.to_spherical that returns an array of [rho, theta, phi] around `center`.
        group_gas_distance_from_center = g3.to_spherical(group_gas_data['POS '], center).T[0]
        group_gas_mask = group_gas_distance_from_center < r200c
        
        
        group_gas_masswtemp =  group_gas_data['TEMP'][group_gas_mask] * group_gas_data['MASS'][group_gas_mask]
        group_gas_avgtemp =  np.mean(group_gas_masswtemp)/np.sum( group_gas_data['MASS'][group_gas_mask])

        print('   Avg Temp(R<R200c): ', group_gas_avgtemp, '[K]')
