# %%
import xtrack as xt
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# add current folder to path
# import sys
# sys.path.append(os.getcwd())


import configure_and_track as configure_and_track

# %%
collider = xt.Multiline.from_json('collider.json')


# %%
collider.build_trackers()

# %%
# collider.vars['beambeam_scale'] = 0 
twiss_b1 = collider['lhcb1'].twiss()
twiss_b2 = collider['lhcb2'].twiss().reverse()

# %%
survey_b1 = {}
survey_b2 = {}

for my_ip in [1,2,5,8]:
    print(f'Survey for IP{my_ip}...')
    survey_b1[f'ip{my_ip}'] = collider['lhcb1'].survey(element0=f'ip{my_ip}')
    survey_b2[f'ip{my_ip}'] = collider['lhcb2'].survey(element0=f'ip{my_ip}').reverse()
# collider.vars['beambeam_scale'] = 1

# %% filling scheme computation
config, config_sim, config_collider = configure_and_track.read_configuration()

filling_scheme = (config_collider['config_beambeam']
                                 ['mask_with_filling_pattern']
                                 ['pattern_fname'])

b1_bunch_to_track = (config_collider['config_beambeam']
                                 ['mask_with_filling_pattern']
                                 ['i_bunch_b1'])
b2_bunch_to_track = (config_collider['config_beambeam']
                                 ['mask_with_filling_pattern']
                                 ['i_bunch_b2'])

# %%

import fillingpatterns as fp
bb_schedule = fp.FillingPattern.from_json(filling_scheme)
bb_schedule.b1.n_bunches
bb_schedule.b2.n_bunches

bb_schedule.n_coll_ATLAS
bb_schedule.n_coll_LHCb
bb_schedule.n_coll_ALICE

bb_schedule.compute_beam_beam_schedule(
    n_lr_per_side=25)    

for ii,zz in zip([bb_schedule.b1,bb_schedule.b2],['Beam 1','Beam 2']):
    my_bb_schedule= ii.bb_schedule.sort_values(by=['collides in ATLAS/CMS',
                                                'collides in LHCB',
                                                'collides in ALICE',
                                                '# of LR in ATLAS/CMS', 
                                                '# of LR in ALICE', 
                                                '# of LR in LHCB',
                                                ], ascending=False)   

    print(f'Suggested bunch ID for {zz}: {my_bb_schedule.index[0]}') 
# %%
bb_schedule_b1 = bb_schedule.b1.bb_schedule.loc[b1_bunch_to_track]
bb_schedule_b2 = bb_schedule.b2.bb_schedule.loc[b2_bunch_to_track]

print('\nBunch to track in Beam 1:')
print(bb_schedule_b1)
print('\nBunch to track in Beam 2:')
print(bb_schedule_b2)

# %% Compute the luminosity
from xtrack import lumi
assert twiss_b1.T_rev0 == twiss_b1.T_rev0

for ii, colliding_bunches in zip(['ip1','ip2','ip5','ip8'],
                                [1,
                                 0,
                                 1,
                                 0]):
    aux = lumi.luminosity_from_twiss(
        colliding_bunches,
        config_collider['config_beambeam']['num_particles_per_bunch'],
        ii,  
        config_collider['config_beambeam']['nemitt_x'],
        config_collider['config_beambeam']['nemitt_y'],
        config_collider['config_beambeam']['sigma_z'],
        twiss_b1,
        twiss_b2,
        crab=False,                          
    )

    sigma_tot = 81e-27 # cm^2
    print(f'Luminosity in {ii}: {aux:.2e} cm^-2 s^-1')
    # compute pile-up from luminosity
    print(f'Pile-up in {ii}: {aux*sigma_tot/colliding_bunches*twiss_b1.T_rev0:.2e}\n')
# %%
for my_ip in ['on_alice_normalized','on_lhcb_normalized']:
    print(f'*****************\nValues for {my_ip} (polarity):')
    print(collider.vars[my_ip]._value)
    print(f'*****************\n')

# %%
collider.vars['beambeam_scale'] = 1
collider.vars['bb_ho.l2b2_05_scale_strength'] = 1 #deactivate (for example) this lens
twiss_b1 = collider['lhcb1'].twiss()
twiss_b2 = collider['lhcb2'].twiss().reverse()
for my_ip in [1,2,5,8]:
    print(f'*****************\nValues for IP{my_ip}:')
    my_df = []
    for my_table, my_beam in zip([twiss_b1, twiss_b2],[1,2]):
        my_df.append(my_table[
            ['x', 'y', 'px', 'py', 'betx', 'bety', 'alfx', 'alfy'],
            f'ip{my_ip}'].to_pandas())
        my_df[-1].index = [f'B{my_beam}']
    print(pd.concat(my_df, axis=0).transpose())
    print(f'*****************\n')

# %%
for ii in collider.vars.get_independent_vars():
    if 'bb_ho' in ii:
        print(ii)
# %%
collider.vars['beambeam_scale']._value 
# %%
for ii in config_collider['config_knobs_and_tuning']['knob_settings'].keys():
    if len(collider.vars[ii]._find_dependant_targets())==1:
        print(ii)
# %%
plt.plot(twiss_b1['s'], twiss_b1['x'], label='x')
plt.plot(twiss_b1['s'], twiss_b1['y'], label='y')
plt.title('Beam 1 closed orbit')
plt.xlabel('s [m]')
plt.ylabel('x, y [m]')
plt.legend()
plt.grid(True)
# %%
def set_orbit_flat(collider):
    print('Setting optics as flat')
    for ii in ['on_x1', 'on_sep1', 'on_x2h', 'on_sep2h', 'on_x2v', 'on_sep2v', 'on_x5', 
               'on_sep5', 'on_x8h', 'on_sep8h', 'on_x8v', 'on_sep8v', 'on_disp', 
               'on_alice_normalized', 'on_lhcb_normalized','on_sol_atlas', 'on_sol_cms', 
               'on_sol_alice', 'i_oct_b1', 'i_oct_b2']:
        collider.vars[ii] = 0
        


set_orbit_flat(collider)
twiss_b1 = collider['lhcb1'].twiss()
twiss_b2 = collider['lhcb2'].twiss().reverse()
plt.plot(twiss_b1['s'], twiss_b1['x'])
plt.plot(twiss_b1['s'], twiss_b1['y'])
plt.title('Beam 1 closed orbit')
plt.xlabel('s [m]')
plt.ylabel('x, y [m]')
plt.grid(True)
# %%
config, config_sim, config_collider = configure_and_track.read_configuration()


def set_orbit_from_config(collider, config):
    print('Setting optics as from config')
    for ii in ['on_x1', 'on_sep1', 'on_x2h', 'on_sep2h', 'on_x2v', 'on_sep2v', 'on_x5', 
               'on_sep5', 'on_x8h', 'on_sep8h', 'on_x8v', 'on_sep8v', 'on_disp', 
               'on_alice_normalized', 'on_lhcb_normalized', 'on_sol_atlas', 'on_sol_cms', 
               'on_sol_alice', 'i_oct_b1', 'i_oct_b2']:
        collider.vars[ii] = config['config_collider']['config_knobs_and_tuning']['knob_settings'][ii]

set_orbit_from_config(collider, config)
twiss_b1 = collider['lhcb1'].twiss()
twiss_b2 = collider['lhcb2'].twiss().reverse()
plt.plot(twiss_b1['s'], twiss_b1['x'])
plt.plot(twiss_b1['s'], twiss_b1['y'])
plt.grid(True)
# %%
#removing bb to understand beta-beating
collider.vars['beambeam_scale'] = 0
twiss_b1 = collider['lhcb1'].twiss()
twiss_b2 = collider['lhcb2'].twiss().reverse()
plt.plot(twiss_b1['s'], np.sqrt(twiss_b1['betx']))
plt.plot(twiss_b1['s'], np.sqrt(twiss_b1['bety']))
plt.grid(True)
#putting beam beam back on
collider.vars['beambeam_scale'] = 1
twiss_b1 = collider['lhcb1'].twiss()
twiss_b2 = collider['lhcb2'].twiss().reverse()
# %%

# with collider.vars['beambeam_scale']:
#     print(collider.vars['beambeam_scale']._value)
#     twiss_b1 = collider['lhcb1'].twiss()
#     twiss_b2 = collider['lhcb2'].twiss().reverse()
#     plt.plot(twiss_b1['s'], np.sqrt(twiss_b1['betx']))
#     plt.plot(twiss_b1['s'], np.sqrt(twiss_b1['bety']))
print(collider.vars['beambeam_scale']._value)

# %%

#twiss_b2[['s','dx_zeta','dy_zeta'],'.*BPLV.*']

set_orbit_from_config(collider, config)
collider.vars['beambeam_scale'] = 0
collider.vars['vrf400'] = 1



twiss_b1 = collider['lhcb1'].twiss()
twiss_b2 = collider['lhcb2'].twiss().reverse()
plt.plot(twiss_b1['s'], twiss_b1['dx_zeta']*1e6, label='Crabbing x-z',alpha=0.5)
plt.plot(twiss_b1['s'], twiss_b1['dy_zeta']*1e6, label='Crabbing y-z',alpha=0.5)

crabbing_xz = twiss_b1[['dx_zeta'],'bplh.7r4.b1']*1e6
crabbing_yz = twiss_b1[['dy_zeta'],'bplv.a6r4.b1']*1e6
plt.axvline(twiss_b1[['s'],'bplh.7r4.b1'], label=f'Crabbing x-z at HT bplh.7r4.b1  ={np.round(crabbing_xz,1)}'+r' $\mu$rad', color='b', linestyle='--')
plt.axvline(twiss_b1[['s'],'bplv.a6r4.b1'], label=f'Crabbing y-z at HT bplv.a6r4.b1  ={np.round(crabbing_yz,1)}'+r' $\mu$rad', color='orange', linestyle='--')
plt.scatter(twiss_b1[['s'],'bplh.7r4.b1'], crabbing_xz, color='b')
plt.scatter(twiss_b1[['s'],'bplv.a6r4.b1'], crabbing_yz, color='orange')
plt.title('Crabbing along the lattice')
plt.legend()
plt.ylabel(r'dx_zeta,dy_zeta [rad]')
plt.xlabel('s [m]')
plt.xlim(twiss_b1[['s'],'bplv.a6r4.b1']-100,twiss_b1[['s'],'bplv.a6r4.b1']+100)
plt.grid(True)

# %%
set_orbit_from_config(collider, config)
collider.vars['beambeam_scale'] = 0
all_crabbing_xz = []
all_crabbing_yz = []
for ii in range(24):
    collider.vars['vrf400'] = ii+1e-3
    twiss_b1 = collider['lhcb1'].twiss()
    crabbing_xz = twiss_b1[['dx_zeta'],'bplh.7r4.b1']*1e6
    crabbing_yz = twiss_b1[['dy_zeta'],'bplv.a6r4.b1']*1e6
    all_crabbing_xz.append(crabbing_xz)
    all_crabbing_yz.append(crabbing_yz)
all_crabbing_xz = np.array(all_crabbing_xz)
all_crabbing_yz = np.array(all_crabbing_yz)

# %%
plt.plot(np.arange(24), all_crabbing_xz,'-o', label='Crabbing x-z at bplh.7r4.b1',alpha=0.5)
plt.plot(np.arange(24), all_crabbing_yz,'-o', label='Crabbing y-z at bplv.a6r4.b1',alpha=0.5)
plt.title('Crabbing as a function of vrf400')
plt.axvline(12, color='red', linestyle='--', label='Nominal vrf400')
plt.grid()
plt.legend()
plt.ylabel(r'dx_zeta,dy_zeta [rad]')
plt.xlabel('vrf400 [MV]')

# %%
set_orbit_from_config(collider, config)
collider.vars['on_x1'] = +170
collider.vars['vrf400'] = 12
all_crabbing_xz = []
all_crabbing_yz = []
beambeam_scale = np.linspace(0,2.5,24)
for ii in range(24):
    collider.vars['beambeam_scale'] = beambeam_scale[ii]
    twiss_b1 = collider['lhcb1'].twiss()
    crabbing_xz = twiss_b1[['dx_zeta'],'bplh.7r4.b1']*1e6
    crabbing_yz = twiss_b1[['dy_zeta'],'bplv.a6r4.b1']*1e6
    all_crabbing_xz.append(crabbing_xz)
    all_crabbing_yz.append(crabbing_yz)
all_crabbing_xz = np.array(all_crabbing_xz)
all_crabbing_yz = np.array(all_crabbing_yz)
# %%
plt.plot(beambeam_scale*1.6, all_crabbing_xz,'-o', label='Crabbing x-z at bplh.7r4.b1',alpha=0.5)
plt.plot(beambeam_scale*1.6, all_crabbing_yz,'-o', label='Crabbing y-z at bplv.a6r4.b1',alpha=0.5)
plt.title('Crabbing as a function of number of particles per bunch')
plt.axvline(1*1.6, color='red', linestyle='--', label='Nominal particle per bunch 1.6e11')
plt.grid()
plt.legend()
plt.ylabel(r'dx_zeta,dy_zeta [rad]')
plt.xlabel('Number of particles per bunch')

# %%
set_orbit_from_config(collider, config)
collider.vars['vrf400'] = 12
all_crabbing_xz = []
all_crabbing_yz = []
crossing_5 = np.linspace(-300,300,24)
for ii in range(24):
    collider.vars['on_x5'] = crossing_5[ii]
    twiss_b1 = collider['lhcb1'].twiss()
    crabbing_xz = twiss_b1[['dx_zeta'],'ip5']*1e6
    crabbing_yz = twiss_b1[['dy_zeta'],'ip5']*1e6
    all_crabbing_xz.append(crabbing_xz)
    all_crabbing_yz.append(crabbing_yz)
all_crabbing_xz = np.array(all_crabbing_xz)
all_crabbing_yz = np.array(all_crabbing_yz)
# %%
plt.plot(crossing_5, all_crabbing_xz,'-o', label='Crabbing x-z at ip5',alpha=0.5)
plt.plot(crossing_5, all_crabbing_yz,'-o', label='Crabbing y-z at ip5',alpha=0.5)
plt.title('Crabbing as a function of on_x5')
plt.axvline(170, color='red', linestyle='--', label='Nominal on_x5')
plt.grid()
plt.legend()
plt.ylabel(r'dx_zeta,dy_zeta [rad]')
plt.xlabel('Number of particles per bunch')

# %%
set_orbit_from_config(collider, config)
collider.vars['vrf400'] = 12
all_crabbing_xz = []
all_crabbing_yz = []
crossing_1 = np.linspace(-300,300,24)
for ii in range(24):
    collider.vars['on_x1'] = crossing_1[ii]
    twiss_b1 = collider['lhcb1'].twiss()
    crabbing_xz = twiss_b1[['dx_zeta'],'ip1']*1e6
    crabbing_yz = twiss_b1[['dy_zeta'],'ip1']*1e6
    all_crabbing_xz.append(crabbing_xz)
    all_crabbing_yz.append(crabbing_yz)
all_crabbing_xz = np.array(all_crabbing_xz)
all_crabbing_yz = np.array(all_crabbing_yz)
# %%
plt.plot(crossing_5, all_crabbing_xz,'-o', label='Crabbing x-z at ip1',alpha=0.5)
plt.plot(crossing_5, all_crabbing_yz,'-o', label='Crabbing y-z at ip1',alpha=0.5)
plt.title('Crabbing as a function of on_x1')
plt.axvline(-170, color='red', linestyle='--', label='Nominal on_x1')
plt.grid()
plt.legend()
plt.ylabel(r'dx_zeta,dy_zeta [rad]')
plt.xlabel('Number of particles per bunch')



# %%

#plt.plot(twiss_b1['s'], twiss_b1['dx'], label='Dispersion x',alpha=0.5)
plt.plot(twiss_b1['s'], twiss_b1['dy'], label='Dispersion y',alpha=0.5)


plt.title('Dispersion along the lattice')
plt.legend()
plt.ylabel(r'dx,dy')
plt.xlabel('s [m]')
plt.grid(True)
# %%
# %%
import xpart as xp
collider.vars['beambeam_scale'] = 0
def set_orbit_flat(collider):
    print('Setting optics as flat')
    for ii in ['on_x1', 'on_sep1', 'on_x2h', 'on_sep2h', 'on_x2v', 'on_sep2v', 'on_x5', 
               'on_sep5', 'on_x8h', 'on_sep8h', 'on_x8v', 'on_sep8v', 'on_disp', 
               'on_alice_normalized', 'on_lhcb_normalized','on_sol_atlas', 'on_sol_cms', 
               'on_sol_alice', 'i_oct_b1', 'i_oct_b2']:
        
        collider.vars[ii] = 0
        collider.vars['on_lhcb_normalized'] = 0.0
        collider.vars['on_x1'] = 0
collider.vars['vrf400'] = 12
set_orbit_flat(collider)
#particle_0 = xp.Particles(mass0=xp.PROTON_MASS_EV, q0=1, p0c=450e9, x=0, px =0, y=0, py=0, zeta=0, delta=0)
#twiss_init = 'periodic', ele_start = 'ip8', ele_stop = 'bb_ho.c8b1_00'
# get_twiss_init(self, at_element)
# twiss_init_b1 = collider['lhcb1'].get_twiss_init()
twiss_b1 = collider['lhcb1'].twiss()
# twiss_b1 = collider['lhcb1'].twiss(twiss_init = twiss_b1.get_twiss_init('ip2'), ele_start = 'drift_12055_part0', ele_stop = 'ip2')
twiss_b2 = collider['lhcb2'].twiss().reverse()

# %%
# start = twiss_b1[:, 'ip8%%-500': 'ip8%%-500']['name'][0]
# end = twiss_b1[:, 'ip8%%+500': 'ip8%%+500']['name'][0]
start = twiss_b2[:, 'ip8%%+0': 'ip8%%+0']['name'][0]
end = twiss_b2[:, 'ip2%%+500': 'ip2%%+500']['name'][0]
print(start, end)
#plt.plot(twiss_b1['s', 'ip8': ], twiss_b1['dx', 'ip8': ],label='x')
plt.plot(twiss_b2['s', : ], twiss_b2['dy', : ],label='y')
plt.axvline(twiss_b2[['s'],'ip8'],color = 'black', linestyle='-.', label='IP8')
plt.axvline(twiss_b2[['s'],'ip1'],color = 'green', linestyle='-.', label='IP1')
plt.axvline(twiss_b2[['s'],'ip2'],color = 'red', linestyle='-.', label='IP2')
plt.title(f'Closed orbit with on_lhcb_normalized = {collider.vars["on_lhcb_normalized"]._value}')
# plt.axvline(twiss_b2['s'][11193])
# plt.ylim(-2e-17,2e-17)
# plt.xlim(twiss_b1[:, 'ip8%%-20': 'ip8%%+20']['s'][0], twiss_b1[:, 'ip8%%-20': 'ip8%%+20']['s'][-1])
#plt.axvline(twiss_b1[['s'],'ip1'],color = 'red', linestyle='--')
plt.legend()
plt.grid(True)





# %%
print(collider.vars['dqx.b1_sq']._value)
print(collider.vars['dqy.b1_sq']._value)
# %%

set_orbit_from_config(collider, config)
collider.vars['beambeam_scale'] = 1
collider.vars['vrf400'] = 12


twiss_b1 = collider['lhcb1'].twiss()
twiss_b2 = collider['lhcb2'].twiss().reverse()
plt.plot(twiss_b1['s'], twiss_b1['dx_zeta']*1e6, label='Crabbing x-z',alpha=0.5)
plt.plot(twiss_b1['s'], twiss_b1['dy_zeta']*1e6, label='Crabbing y-z',alpha=0.5)

crabbing_xz = twiss_b1[['dx_zeta'],'bplh.7r4.b1']*1e6
crabbing_yz = twiss_b1[['dy_zeta'],'bplv.a6r4.b1']*1e6
plt.axvline(twiss_b1[['s'],'bplh.7r4.b1'], label=f'Crabbing x-z at HT bplh.7r4.b1  ={np.round(crabbing_xz,1)}'+r' $\mu$rad', color='r', linestyle='--')
plt.axvline(twiss_b1[['s'],'bplv.a6r4.b1'], label=f'Crabbing y-z at HT bplv.a6r4.b1  ={np.round(crabbing_yz,1)}'+r' $\mu$rad', color='r', linestyle='--')
plt.title('Crabbing along the lattice and at the HT monitors')
plt.legend()
plt.ylabel(r'dx_zeta,dy_zeta [$\mu$ rad]')
plt.xlabel('s [m]')
plt.grid(True)
#get the s at 
# %%
collider_start = xt.Multiline.from_json('/afs/cern.ch/work/a/afornara/public/run3/example_DA_study/master_study/master_jobs/1_build_distr_and_collider/collider/collider.json')

# %%
collider_start.build_trackers()

# %%
set_orbit_flat(collider)
# collider.vars['beambeam_scale'] = 0 
twiss_b1_start = collider_start['lhcb1'].twiss()
twiss_b2_start = collider_start['lhcb2'].twiss().reverse()



# %%
