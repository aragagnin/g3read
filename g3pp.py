"""

This file compute gravitational potential on a simcut snapshot.

Edit the snapshot folder and properties at the end of this file.

You can test this library using the function test()

You can also use this file as 'library' and call gravitational_potential

Antonio Ragagnin

"""

import g3read as g3r, numpy as np, sys, math, scipy.stats

def gravitational_potential(masses, #array of masses
                            positions = None, # array of particle x,y,z positions or r,theta,phi
                            spherical = False, # if true: read positions as spherical coordinates
                            gpos = None, # if spherical is False then you must specify the center of the halo used to compute relative position
                                         # (I suppose positions are NOT relative to the halo)
                            cut = None, #if not None, compute potential only within a radius of value `cut`
                            G=43007.1, #gravitational G constant, default is in Gadget code unites
                            spher_nbs=40, #number of binning in the spherical radius coordinate
                            spher_nfi=4,  #number of binning in the spherical phi coordinate
                            spher_nteta=4  #number of binning in the spherical theta coordinate
                            ):

    all_data={}
    all_data["MASS"]=masses

    if spherical == False:
        all_sferical = g3r.to_spherical(positions, gpos)
    else:
        all_sferical = positions

    all_data["SPOS"] = all_sferical

    Nall=len(all_data['MASS'])



    twopi=2.*math.pi
    pi=math.pi

    """    POTENTIAL    """

    spher_bs = [np.logspace(np.log10(np.min(all_data['SPOS'][:,0])+0.01),np.log10(np.max(all_data['SPOS'][:,0])),spher_nbs),np.linspace(0.,pi,spher_nteta), np.linspace(-pi,pi,spher_nfi)]


    mass_weights = all_data['MASS']
    if cut is not None:
        mass_weights[all_data['SPOS'][:,0]>cut]=0.

    spher_all_ms, spher_b = np.histogramdd(all_data['SPOS'], weights=mass_weights, bins=spher_bs)

    spher_all_ds, spher_b = np.histogramdd(all_data['SPOS'], weights=all_data['SPOS'].T[0], bins=spher_bs)
    spher_all_ts, spher_b = np.histogramdd(all_data['SPOS'], weights=all_data['SPOS'].T[1], bins=spher_bs)
    spher_all_fs, spher_b = np.histogramdd(all_data['SPOS'], weights=all_data['SPOS'].T[2], bins=spher_bs)
    spher_all_ns, spher_b = np.histogramdd(all_data['SPOS'],  bins=spher_bs)


    spher_all_ns[spher_all_ns==0]=np.nan
    spher_all_cds = spher_all_ds/spher_all_ns
    spher_all_cts = spher_all_ts/spher_all_ns
    spher_all_cfs = spher_all_fs/spher_all_ns


    spher_all_x ,    spher_all_y ,    spher_all_z = g3r.to_cartesian(np.array([spher_all_cds,spher_all_cts,spher_all_cfs]).T)

    shape=spher_all_ds.shape
    spher_b_delta_r=(spher_b[0][1:]-spher_b[0][:-1])
    spher_b_delta_t=(spher_b[1][1:]-spher_b[1][:-1])
    spher_b_delta_f=(spher_b[2][1:]-spher_b[2][:-1])

    shper_delta_rs = np.transpose( (np.transpose(np.ones(shape),axes=(2,1,0) )* (spher_b_delta_r)), axes=(2,1,0))
    shper_delta_ts = np.transpose( (np.transpose(np.ones(shape),axes=(0,2,1) )* (spher_b_delta_t)), axes=(0,2,1))
    shper_delta_fs = np.transpose( (np.transpose(np.ones(shape),axes=(0,1,2) )* (spher_b_delta_f)), axes=(0,1,2))

    spher_all_vols = spher_all_cds**2.*np.sin(spher_all_cts)*shper_delta_rs*shper_delta_ts*shper_delta_fs
    spher_all_rhos = spher_all_ms/spher_all_vols
    spher_all_ms = np.nan_to_num(spher_all_ms)

    def generate_fi(spher_b,spher_all_cds,spher_all_cts,spher_all_cfs,spher_all_x,spher_all_y,spher_all_z,spher_all_ms):
        fi=np.ones(spher_all_ds.shape)
        for bin_r in range(len(spher_b[0])-1):
            for bin_t in range(len(spher_b[1])-1):
                for bin_phi in range(len(spher_b[2])-1):
                    position_xyz = g3r.to_cartesian(np.array(np.array([spher_all_cds[bin_r,bin_t,bin_phi], spher_all_cts[bin_r,bin_t,bin_phi],spher_all_cfs[bin_r,bin_t,bin_phi]])).T)
                    distances = np.sqrt( (spher_all_x-position_xyz[0])**2. + (spher_all_y-position_xyz[1])**2. + (spher_all_z-position_xyz[2])**2.)
                    distances = np.nan_to_num(distances)
                    non_zero_distances = distances>0.
                    fi[bin_r,bin_t,bin_phi] = -G*np.sum(spher_all_ms[non_zero_distances]/distances[non_zero_distances])
        return np.nan_to_num(fi)
    fi =  generate_fi(spher_b,spher_all_cds,spher_all_cts,spher_all_cfs,spher_all_x,spher_all_y,spher_all_z,spher_all_ms)


    bin_all_h_i = np.digitize(all_data['SPOS'][:,0],spher_bs[0])-1
    bin_all_h_j = np.digitize(all_data['SPOS'][:,1],spher_bs[1])-1
    bin_all_h_k = np.digitize(all_data['SPOS'][:,2],spher_bs[2])-1

    # bug or weird behviour of  of np:
    # if a value is exactly a boundary, the bin is larger than it should
    bin_all_h_i[ bin_all_h_i>=len(spher_bs[0])-1 ]=len(spher_bs[0])-2
    bin_all_h_j[ bin_all_h_j>=len(spher_bs[1])-1 ]=len(spher_bs[1])-2
    bin_all_h_k[ bin_all_h_k>=len(spher_bs[2])-1 ]=len(spher_bs[2])-2

    bin_all_h=np.array([bin_all_h_i,bin_all_h_j,bin_all_h_k]).T

    somma_all_inte = fi[tuple ( bin_all_h.T)]


    return somma_all_inte



def test():
    my_file="nitk/new_snap_136"
    center_x = 208791.203
    center_y = 203941.000
    center_z = 349847.000
    r500_cu= 2517.0
    radial_bins=30


    center = np.array([center_x, center_y, center_z])
    data = g3r.read_new(my_file, ['MASS',  'POS '], -1)

    #it convert an array of [x,y,z] -> [r, theta, phi]
    data['SPOS'] = g3r.to_spherical(data['POS '], center)

    pot = gravitational_potential(data['MASS'],
                                  positions = data['SPOS'],
                                  spherical = True,
                                  cut = r500_cu
    )


    h_s, h_b, h_i = scipy.stats.binned_statistic(data['SPOS'][:,0], #component zero is radius
                                                 pot,
                                                 bins = np.logspace(np.log10(10.), np.log10(r500_cu), radial_bins),
                                                 statistic = 'mean'
    )

    print('r_min     r_max      pot')
    for i, s in enumerate(h_s):
        print('%.3e %.3e %.3e'%(h_b[i], h_b[i+1] ,s))
