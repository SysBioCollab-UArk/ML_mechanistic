from pysb.examples.earm_1_0 import model
from pysb.simulator import ScipyOdeSimulator
import numpy as np
import matplotlib.pyplot as plt


from scipy.optimize import curve_fit
import os
import read_data
import json
from copy import deepcopy
from read_data import get_k_div

print(model.observables)

print(model.parameters)
#quit()






tspan = np.linspace(0,72*3600,1001)


sim=ScipyOdeSimulator(model,tspan,verbose=True)
#for l_0 in [0.1,0.5,1]:

   # result= sim.run(param_values={'L_0': l_0})


  #  plt.plot(tspan/3600,result.observables['CPARP_total'],lw=2,label='L_0 = %g' % l_0)

#plt.xlabel('time (hr)')
#plt.ylabel('molecule count')
#plt.legend(loc=0)



#plt.show()




#quit()



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


Cell_div_times=list(read_data.get_k_div("data/doubling_times.csv").values())
Cell_key=list(read_data.get_k_div("data/doubling_times.csv").keys())
#Cell_div_times =  [x*4  for x in Cell_div_times]
print(np.array(Cell_key))
print(np.array(Cell_div_times))
#quit()

CellAA = [0]*len(Cell_div_times)
CellIC = [0]*len(Cell_div_times)
vCellAA = [0]*len(Cell_div_times)
vCellIC = [0]*len(Cell_div_times)
ts = np.linspace(0, 72*60*60, 1001)
cell = 0
Progress=0

# loop over all cell lines
for cell_line in Cell_key:   #["HeLa"]:    Cell_key:   ["AN3-CA"]:

    # print(t_div[cell_line])
    Progress = Progress + 1
    print(Progress)
    # set initial concentrations for each protein

    gene_expr_ratio = read_data.get_gene_expr(os.path.join('data', 'normalized_ratios.csv'), cell_line)
    initials = {}
    for gene in gene_map.keys():
        for p_name in gene_map[gene]:
            idx = np.where(np.array([ic[1].name for ic in model.initial_conditions]) == p_name)[0][0]
            sp = model.initial_conditions[idx][0]
            amount = round(model.parameters[p_name].value * gene_expr_ratio[gene])
            initials[sp] = amount
            print(amount)
            #quit()

    #print(Cell_div_times)
    #Cell_div_times=np.array([0.2,0.13524823035316005,0.0009973340727481227, 0.01933464938800405, 0.015322684707684294, 0.009560650766344074,0.024755256448569473,0.017503716680806698,0.0204669443078724,0.016681350442654508,0.021061901566695386,0.02161132344294155,0.016675553662228032 ,0.013119504994825462,0.011375500775053806,0.011108127893588867,0.015576341136178546,0.018216745875425634])


    #print(Cell_div_times)

    #quit()
    sim = ScipyOdeSimulator(model, ts)


    #print (model.observables)


    n = 10                             #30

    Kprolif = [0]*n
    TTD = [0]*n
    ProbL = [0]*n
    Start = 10**np.linspace(1, 8, n)         #10**np.linspace(-15, 8, n)
    K_Div = 40      #Cell_div_times[cell]/3600    # time units 1/seconds
    #quit()


    Count = 0
    print(Start)



    while Count < n:
        #newinitial = deepcopy(initials)
        print("Count",Count)
        idx = np.where(np.array([ic[1].name for ic in model.initial_conditions]) == "L_0")[0][0]
        sp = model.initial_conditions[idx][0]


        initials[sp] = Start[Count] #if Count==0 else np.array([Start[Count]])

        #initials[model.monomers["L"](b=None)] = Start[Count] if Count == 0 else np.array([Start[Count]])
        print("start count",Start[Count])
        print(initials)

        traj = sim.run(initials= initials)
        traj = sim.run(initials=initials, param_values={'kc8': 0})
        TODT = 0

        #print(Start)
        #quit()
        ######################################################################################################
        #print(traj.dataframe["CPARP_total"])

        #plt.plot(traj.dataframe["CPARP_total"])

        #plt.show()
       # quit()
        P = 0
        Q = 0
        print(max(traj.dataframe["CPARP_total"]))
        for T in traj.dataframe["CPARP_total"]:
            # print(T)
            # print(P)
            P = P + 1
            if T < model.parameters["PARP_0"].value/2:
                TODF = T
                D = P - 1
            if T > model.parameters["PARP_0"].value/2 and Q == 0:
                TODT = T
                Q = Q + 1

        # print("TODF", TODF)
        # print("TODT", TODT)
        # print("D", D)
        #TOD = D * 200                             #* 200
        #print("TOD", TOD)

        Ya = D
        Yb = D + 1
        Xa = TODF
        Xb = TODT
        X = model.parameters["PARP_0"].value/2
        Y = 1
        Y = Ya + (X - Xa) * ((Yb - Ya) / (Xb - Xa))
        print("Y =", Y)
        TOD = Y * 259.2

        ########################################################################################################

        TTD[Count] = TOD
        Death = TOD     # Time of death
        domain = ts        # timeScale
        x = domain



        print('x=',x)
        print('Death=',Death)
        #quit()
        ##############################################################################################################






        plt.plot(ts / 3600, traj.observables["CPARP_total"], color="0.5")
        plt.title(Cell_key[cell])
        plt.xlabel('time (h)')
        plt.ylabel('concentration')
        plt.axvline(x = Death/3600)
        #plt.show()

        Count = Count+1
    plt.show()
    cell = cell + 1




plt.show()

