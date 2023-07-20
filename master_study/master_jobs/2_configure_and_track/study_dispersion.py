# %%
import xtrack as xt
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# %%
collider = xt.Multiline.from_json('/afs/cern.ch/work/a/afornara/public/run3/example_DA_study'
                                  '/master_study/master_jobs/2_configure_and_track/collider.json')

# %%
collider.build_trackers()

# %%

def set_orbit_flat(collider):
    print('Setting optics as flat')
    for ii in ['on_x1', 'on_sep1', 'on_x2h', 'on_sep2h', 'on_x2v', 'on_sep2v', 'on_x5', 
               'on_sep5', 'on_x8h', 'on_sep8h', 'on_x8v', 'on_sep8v', 'on_disp', 
               'on_alice_normalized', 'on_lhcb_normalized','on_sol_atlas', 'on_sol_cms', 
               'on_sol_alice', 'i_oct_b1', 'i_oct_b2']:
        collider.vars[ii] = 0

def set_orbit_from_config(collider, config):
    print('Setting optics as from config')
    for ii in ['on_x1', 'on_sep1', 'on_x2h', 'on_sep2h', 'on_x2v', 'on_sep2v', 'on_x5', 
               'on_sep5', 'on_x8h', 'on_sep8h', 'on_x8v', 'on_sep8v', 'on_disp', 
               'on_alice_normalized', 'on_lhcb_normalized', 'on_sol_atlas', 'on_sol_cms', 
               'on_sol_alice', 'i_oct_b1', 'i_oct_b2']:
        collider.vars[ii] = config['config_collider']['config_knobs_and_tuning']['knob_settings'][ii]
# %%
#First problem: the x closed orbit starts at zero and then an oscillation begins at IP8 and
#the y closed orbit is oscillating until IP2
collider.vars['beambeam_scale'] = 0
collider.vars['vrf400'] = 12
set_orbit_flat(collider)
twiss_b1 = collider['lhcb1'].twiss()
twiss_b2 = collider['lhcb2'].twiss().reverse()


# %% It seems that the HO is kicking...
# It takes a while to run (2-3 minutes)
import xpart as xp
my_list = []
for ii in collider['lhcb1'].element_dict:
    ##if 'bb' in ii:
        my_particle = xp.Particles(p0c=450e9, #eV
             q0=1,
             mass0=xp.PROTON_MASS_EV)
        collider['lhcb1'].element_dict[ii].track(my_particle)
        my_list.append({'name':ii,'delta_px':my_particle.px[0],'delta_py':my_particle.py[0]})
        #print(f'element {ii}: delta_px = {my_particle.px[0]},  delta_py = {my_particle.py[0]}')
my_df = pd.DataFrame(my_list)
my_df[my_df['delta_px'] != 0]
my_df[my_df['delta_py'] != 0]
# %%
import xpart as xp
my_particle = xp.Particles(p0c=450e9, #eV
             q0=1, delta=1e-4,
             mass0=xp.PROTON_MASS_EV)
print(my_particle.to_dict())
collider['lhcb1'].element_dict['mbxwt.1l2'].track(my_particle)
print(my_particle.to_dict())
# %% 
collider['lhcb1'].element_dict['bb_lr.r1b1_01'].to_dict()
# %%
collider['lhcb1'].element_dict['bb_ho.c1b1_00'].to_dict()

# %%
plt.plot(twiss_b1['s', : ], twiss_b1['x', : ],label='x')
plt.plot(twiss_b2['s', : ], twiss_b2['y', : ],label='y')
plt.axvline(twiss_b2[['s'],'ip8'],color = 'green', linestyle='-.', label='IP8')
plt.axvline(twiss_b2[['s'],'ip1'],color = 'black', linestyle='-.', label='IP1')
plt.axvline(twiss_b2[['s'],'ip2'],color = 'red', linestyle='-.', label='IP2')
plt.title(f'Closed orbit')
plt.ylim(-2e-17,2e-17)
plt.legend()
plt.grid(True)
# %%
# %%
#Second problem: y dispersion
collider.vars['beambeam_scale'] = 0
collider.vars['vrf400'] = 12
set_orbit_flat(collider)
twiss_b1 = collider['lhcb1'].twiss()
twiss_b2 = collider['lhcb2'].twiss().reverse()

# %%
#plt.plot(twiss_b1['s', : ], twiss_b1['dx', : ],label='x')
plt.plot(twiss_b2['s', : ], twiss_b2['dy', : ],label='y')
plt.axvline(twiss_b2[['s'],'ip8'],color = 'green', linestyle='-.', label='IP8',alpha = 0.5)
plt.axvline(twiss_b2[['s'],'ip1'],color = 'black', linestyle='-.', label='IP1',alpha = 0.5)
plt.axvline(twiss_b2[['s'],'ip2'],color = 'red', linestyle='-.', label='IP2',alpha = 0.5)
plt.title(f'Dispersion')
plt.legend()
plt.grid(True)
# %%
#We can study the rotated lattice
collider_start = xt.Multiline.from_json('/afs/cern.ch/work/a/afornara/public/run3/example_DA_study/master_study/master_jobs/1_build_distr_and_collider/collider/collider.json')

# %%
collider_start.build_trackers()

# %%
#Here the orbit is zero!
set_orbit_flat(collider_start)
# collider.vars['beambeam_scale'] = 0 
twiss_b1_start = collider_start['lhcb1'].twiss()
twiss_b2_start = collider_start['lhcb2'].twiss().reverse()

plt.plot(twiss_b1_start['s', : ], twiss_b1_start['x', : ],label='x')
plt.plot(twiss_b2_start['s', : ], twiss_b2_start['y', : ],label='y')
plt.axvline(twiss_b2_start[['s'],'ip8'],color = 'green', linestyle='-.', label='IP8')
plt.axvline(twiss_b2_start[['s'],'ip1'],color = 'black', linestyle='-.', label='IP1')
plt.axvline(twiss_b2_start[['s'],'ip2'],color = 'red', linestyle='-.', label='IP2')
plt.title(f'Closed orbit')
plt.ylim(-2e-17,2e-17)
plt.legend()
plt.grid(True)
# %%

set_orbit_flat(collider_start)
# collider.vars['beambeam_scale'] = 0 
twiss_b1_start = collider_start['lhcb1'].twiss()
twiss_b2_start = collider_start['lhcb2'].twiss().reverse()

#plt.plot(twiss_b1_start['s', : ], twiss_b1_start['dx', : ],label='x')
plt.plot(twiss_b2_start['s', : ], twiss_b2_start['dy', : ],label='y')
plt.axvline(twiss_b2_start[['s'],'ip8'],color = 'green', linestyle='-.', label='IP8')
plt.axvline(twiss_b2_start[['s'],'ip1'],color = 'black', linestyle='-.', label='IP1')
plt.axvline(twiss_b2_start[['s'],'ip2'],color = 'red', linestyle='-.', label='IP2')
plt.title(f'Dispersion')
plt.legend()
plt.grid(True)
# %%
aux = twiss_b1[['dy','name'],:].to_pandas()
aux[aux['dy']>1e-8]
# %%
aux = collider['lhcb1'].twiss(ele_start= 'mbxwt.1l2', 
                              ele_stop='mbxwt.1l2', 
                              twiss_init= 'periodic')
aux['dy',:]
# %%
my_line = xt.Line(
    elements=[xt.Drift(length=2.), collider['lhcb1'].element_dict['mbxwt.1l2']],
    element_names=['drift_0', 'x'])
my_line.particle_ref = xp.Particles(
                    mass0=xp.PROTON_MASS_EV, q0=1, energy0=0.450e12)
my_line.build_tracker()
my_line.twiss(twiss_init='preserve')
# %%
# %%
import xpart as xp
my_particle = xp.Particles(p0c=450e9, #eV
             q0=1, delta=1e-4,
             mass0=xp.PROTON_MASS_EV)
print(my_particle.to_dict())
collider['lhcb1'].track(my_particle)

# %%
# %%
collider['lhcb1'].discard_tracker()

my_index = collider['lhcb1'].element_names.index('mbxwt.1l2')
my_index = collider['lhcb1'].element_names.index('ip5')


my_line = xt.Line(
    elements=collider['lhcb1'].element_dict,
    element_names=[ii for ii in collider['lhcb1'].element_names][:my_index-1])
my_line.particle_ref = xp.Particles(
                    mass0=xp.PROTON_MASS_EV, q0=1, energy0=0.450e12)
my_line.build_tracker()
aux = xp.Particles(
                    mass0=xp.PROTON_MASS_EV, q0=1, energy0=0.450e12, delta=1e-4)
my_line.track(aux)
print(aux.to_dict()['y'])
# %%
# get position of ip3 in the list of elements
ip3_index = collider['lhcb1'].element_names.index('ip3')
