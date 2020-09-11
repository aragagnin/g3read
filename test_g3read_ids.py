"""

This file tests g3matcha library to loop over haloes, sub haloes and particle extractions combined with g3read.

Note: `snapbase` and `groupbase` values are meant to work on dorc1.usm.uni-muenchen.de servers.
change it by hand to test it in your environment.

Antonio Ragagnin

"""

import g3read as g3, g3matcha as matcha, numpy as np

""" INPUT TEST PARAMETER TO READ MAGNETICUM SIMS """

groupbase = '/HydroSims/Magneticum/Box1a/mr_bao/groups_144/sub_144'

#if you do not need IDs, set it to false for a faster lookup
with_ids = True

#if you do not need to obtain particles of each halo/subhalo,
#set to False for a faster execution time
read_particles = True
snapbase = '/HydroSims/Magneticum/Box1a/mr_bao/snapdir_144/snap_144'

#if you read particles of haloes and subhaloes, set here the quantities you are interested in
read_particles_blocks = ['MASS','ID  ','POS ']

#switch it on to compute concentration of halo 
compute_concentration = True

""" BEGIN """

printf = g3.printf

subhalo_ifile=0

#
# Here we loop through haloes in the catalog.
# Note that blocks must be given in tuples (i.e. with round brakets).
#
for halo  in  matcha.yield_haloes(groupbase, with_ids=with_ids, blocks=('GLEN', 'MCRI', 'RCRI', 'GPOS')):
    ihalo = halo['ihalo']
    #
    # we display properties of the current haloes
    #
    printf('- ihalo: %d\n'% ihalo)
    printf('  GLEN: %d\n'% halo['GLEN'])
    printf('  MCRI [1e10 Msun]: %d\n'% halo['MCRI'])
    printf('  RCRI [akpc/h]: %d\n'% halo['RCRI'])
    printf('  GPOS [akpc/h]: [%.1f, %.1f, %.1f]\n'%tuple(halo['GPOS'].tolist()))
    isubhalo =-1

    if with_ids:
        halo_ids = halo['ids']
    else:
        halo_ids = None
    #
    # If required by read_particles, we read all halo particles within RCRI with read_particles_in_box
    #
    if read_particles:
        halo_particles = g3.read_particles_in_box(snapbase, halo['GPOS'], halo['RCRI'], read_particles_blocks, -1)
        printf('  Readed n. of halo particles: %d\n'%len(halo_particles['MASS']))
        if compute_concentration:
            halo_particles['DIST'] = g3.to_spherical(halo_particles['POS '], halo['GPOS']).T[0] #get distnaces from center
            nfw_res = matcha.nfw_fit_fast_cu(halo_particles['MASS'], halo_particles['DIST'], halo['RCRI'])
            printf('  c200c tot. matter : %.1f\n'% nfw_res['c'])
            
    printf('\n')
    printf('  subhaloes: \n')
    #
    # Here we loop over subhaloes of the current halo
    # Note we must pass halo['GOFF'] to the subhalo-reader in order to efficiently find subhaloes.
    #
    for subhalo in  matcha.yield_subhaloes(groupbase, ihalo, with_ids = with_ids, halo_ids = halo_ids, blocks=('SPOS','MBID','RHMS','MSUB','SMST'), halo_goff = halo['GOFF']):
        isubhalo+=1

        #
        # we only show galaxies within Rcri, if you remove this constrain, then fix the read_particles_in_box call
        #
        if g3.periodic_distance(halo['GPOS'], subhalo['SPOS'], halo['boxsize'])> (halo['RCRI'] - subhalo['RHMS']):
            continue
        
        #
        # we report properties of the current subhalo
        #
        printf('    - isubhalo: %d\n'% isubhalo)
        printf('      SPOS [akpc/h]: [%.1f, %.1f, %.1f]\n'%tuple(subhalo['SPOS'].tolist()))
        printf('      DM half mass radius [akpc/h]: %.1f\n'% subhalo['RHMS'])
        printf('      Total subhalo mass [1e10 Msun/h]: %.1f\n'% subhalo['MSUB'])
        printf('      Stellar mass [1e10 Msun/h]: %.1f\n'% subhalo['SMST'][4])

        subhalo_ifile = subhalo['ifile']
        
        #
        # here we extract sub find particles based on their IDs 
        #
        if read_particles and with_ids:
            # np.in1d finds which of the particles read above are of this subhaloes, based on ID matching
            subhalo_particles_mask = np.in1d( halo_particles['ID  '], subhalo['ids'])
            subhalo_paritcles_mass = halo_particles['MASS'][subhalo_particles_mask]
            printf('      Readed n. of subhalo particles: %d\n'%len(subhalo_paritcles_mass))
            printf('      Total mass from snaphsot particles [1e10 Msun/h]: %.1f\n'% np.sum(subhalo_paritcles_mass))
            

    printf('\n')
