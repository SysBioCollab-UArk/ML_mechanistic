from pysb.examples.earm_1_0 import model
from pysb.simulator import ScipyOdeSimulator
from pysb.simulator import BngSimulator
from read_data import *
import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import splrep, sproot


# get time-to-death
def get_ttd(tspan, cparp_traj):
    cparp_traj_norm = cparp_traj / np.nanmax(cparp_traj)
    st, sc, sk = splrep(tspan, cparp_traj_norm)
    try:
        t10 = sproot((st, sc - 0.10, sk))[0]
        t90 = sproot((st, sc - 0.90, sk))[0]
        t50 = sproot((st, sc - 0.50, sk))[0]
    except IndexError:
        t10 = 0
        t90 = 0
        t50 = 0
    # time-to-death  = halfway point between 10 and 90%
    td = (t10 + t90) / 2
    return (td, t50, cparp_traj_norm[-1])


# dose-response curve
def drc(x, emax, ec, hh):
    return emax + (1 - emax)/(1+(10**x/ec)**hh)


# create the simulator
tspan = np.linspace(0, 20000, 101)
sim = ScipyOdeSimulator(model, tspan, verbose=True)
# sim = BngSimulator(model, tspan, verbose=True)

# for ic in model.initial_conditions:
#     print(ic)
# quit()

# dictionary to map gene names in initial concentrations
gene_map = {
    'BAX': ['Bax_0'],
    'BCL2': ['Bcl2_0', 'Bcl2c_0'],
    'BID': ['Bid_0'],
    'CYCS': ['mCytoC_0'],
    'FADD': ['Bcl2_0'],
    'CASP3': ['pC3_0'],
    'CASP6': ['pC6_0'],
    'CASP8': ['pC8_0'],
    'CASP9': ['pC9_0'],
    'CFLAR': ['flip_0'],  # FLIP
    'XIAP': ['XIAP_0'],
    'DIABLO': ['mSmac_0'],  # SMAC
    'TNFRSF10A': [],  # DR4
    'TNFRSF10B': ['pR_0'],  # DR5
    'PARP1': ['PARP_0'],
    'APAF1': ['Apaf_0'],
    'BFAR': ['BAR_0']
}

# get division times for each cell line
k_div = get_k_div(os.path.join('data', 'doubling_times.csv'))
for kd in k_div.items():
    print(kd)
quit()

# loop over all cell lines
for cell_line in ['HeLa']:  # t_div.keys():
    print(cell_line)
    # print(t_div[cell_line])

    # set initial concentrations for each protein
    gene_expr_ratio = get_gene_expr(os.path.join('data', 'normalized_ratios.csv'), cell_line)
    initials = {}
    for gene in gene_map.keys():
        for p_name in gene_map[gene]:
            idx = np.where(np.array([ic[1].name for ic in model.initial_conditions]) == p_name)[0][0]
            sp = model.initial_conditions[idx][0]
            amount = round(model.parameters[p_name].value * gene_expr_ratio[gene])
            initials[sp] = amount

    # set initial ligand amount
    idx = np.where(np.array([ic[1].name for ic in model.initial_conditions]) == 'L_0')[0][0]
    sp = model.initial_conditions[idx][0]
    initials[sp] = 3000  # 3000

    # run the simulation and calculate time-to-death
    result = sim.run(initials=initials)
    ttd = get_ttd(tspan, result.observables['CPARP_total'])
    print(ttd)

    # plot time courses
    plt.figure()
    plt.title(cell_line)
    for obs in model.observables:
        plt.plot(tspan/3600., result.observables[obs.name], lw=2, label=obs.name)
    plt.xlabel('time (h)')
    plt.ylabel('amount (# molecules)')
    plt.legend(loc=0)

    plt.savefig(os.path.join('plots', '%s.png' % cell_line), format='png')

plt.show()
