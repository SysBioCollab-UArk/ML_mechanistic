from pysb.examples.earm_1_0 import model
from pysb.simulator import ScipyOdeSimulator
from read_data import *
import os
import numpy as np
import matplotlib.pyplot as plt

# create the simulator
tspan = np.linspace(0, 20000, 101)
sim = ScipyOdeSimulator(model, tspan, verbose=True)

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
t_div = get_t_div(os.path.join('data', 'doubling_times.csv'))

# loop over all cell lines
for cell_line in t_div.keys():
    print(cell_line)

    # set initial concentrations for each protein
    gene_expr_ratio = get_gene_expr(os.path.join('data', 'normalized_ratios.csv'), cell_line)
    initials = {}
    for gene in gene_map.keys():
        for p_name in gene_map[gene]:
            idx = np.where(np.array([ic[1].name for ic in model.initial_conditions]) == p_name)[0][0]
            sp = model.initial_conditions[idx][0]
            amount = round(model.parameters[p_name].value * gene_expr_ratio[gene])
            initials[sp] = amount

    # run the simulation and plot time courses
    result = sim.run(initials=initials)
    plt.figure()
    plt.title(cell_line)
    for obs in model.observables:
        plt.plot(tspan, result.observables[obs.name], lw=2, label=obs.name)
    plt.legend(loc=0)

    plt.savefig(os.path.join('plots', '%s.png' % cell_line), format='png')

plt.show()
