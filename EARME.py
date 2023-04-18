from scipy.stats import norm
from scipy.special import erf

from pysb.examples.earm_1_0 import model
#from earm.lopez_embedded import model
import numpy as np
import matplotlib.pyplot as plt
from pysb.simulator import ScipyOdeSimulator
from pysb import *
import math
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
Cell_div_times =  [x*4  for x in Cell_div_times]
print(np.array(Cell_key))
print(np.array(Cell_div_times))
#quit()

CellAA = [0]*len(Cell_div_times)
CellIC = [0]*len(Cell_div_times)
vCellAA = [0]*len(Cell_div_times)
vCellIC = [0]*len(Cell_div_times)
ts = np.linspace(0, 72*60*60, 1001)
cell = 0

# loop over all cell lines
for cell_line in ["HeLa"]:   #["HeLa"]:    Cell_key:   ["AN3-CA"]:

    # print(t_div[cell_line])

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


    n = 30                              #30

    Kprolif = [0]*n
    TTD = [0]*n
    ProbL = [0]*n
    Start = 10**np.linspace(1, 8, n)         #10**np.linspace(-15, 8, n)
    K_Div = Cell_div_times[cell]/3600    # time units 1/seconds
    #quit()


    Count = 0
    print(Start)

    fig, ax1 = plt.subplots()
    ax2 = ax1.twinx()


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
        TOD = Y * 259.2 * 5

        ########################################################################################################

        TTD[Count] = TOD
        Death = TOD     # Time of death
        domain = ts        # timeScale
        x = domain

        ##############################################################################################################

        Prob = 1-np.exp(-K_Div*Death)




        ax1.plot(ts / 3600, traj.observables["CPARP_total"], color="0.5")
        ax1.set_xlabel('time (h)')
        ax1.set_ylabel('concentration')
        if Count == 0 :
            times = []
            probs = []
            for i in range(len(ts) - 1):
                times.append(0.5 * (ts[i] + ts[i + 1]))
                probs.append(1 - np.exp(-K_Div * times[-1]))
            ax2.plot(np.array(times) / 3600, probs, lw=2)

            ax2.set_ylabel('probability')
        plt.tight_layout()











        if TODT == 0:            # if Cparp does not reach 50% there is 100% prob of life
            print("WARNING: CPARP did not reach 50% of maximum")
            Kprolif[Count] = K_Div
            #quit()
        else:

            ProbL[Count] = Prob

            print("Prob of Life =", format(Prob * 100, ".2f"), "%")
            Kprolif[Count] = K_Div * ProbL[Count] - (1. / TTD[Count]) * (1. - ProbL[Count])

        Count += 1

    print("ProbL of Life =", ProbL)
    print("Kprolif =", np.array(Kprolif))
    print("Dip =", np.array(Kprolif)/np.log(2)*3600)
    #for k in Kprolif:
      #  plt.plot(ts / 3600, np.log2(np.exp(ts * k)), label="Kprolif = %g" % k)

    plt.show()
    quit()



    if cell % 3 == 0:
        plt.figure()
        plt.savefig("Viability_Dip_%d.pdf" % cell, format="pdf")

   # DipRate###
    plt.subplot(3, 3, cell%3*3+1)
    plt.title(cell_line)
    plt.plot(ts/3600, np.exp(ts * K_Div), label="Kprolif = %g" % K_Div)
    for k in Kprolif:
        plt.plot(ts/3600, np.exp(ts * k), label="Kprolif = %g" % k)
    plt.xlabel('time (hr)')
    plt.ylabel('cell count')

    #plt.legend(loc="best")

    plt.subplot(3, 3, cell%3*3+2)
    #plt.title("log_tc")
    plt.plot(ts/3600,np.log2(np.exp(ts * K_Div)), label="Kprolif = %g" % K_Div)
    for k in Kprolif:
        plt.plot(ts/3600, np.log2(np.exp(ts*k)), label="Kprolif = %g" % k)
    plt.xlabel('time (hr)')
    plt.ylabel(r'log$_2$ cell count')
    #plt.legend(loc="best")

    #plt.yscale("log", base=2)

    plt.subplot(3,3, cell%3*3+3)
    Dip = np.array(Kprolif)/np.log(2)*3600      # units of 1/hour
    LogConc = np.log10(np.array(Start))



    Dip0=K_Div*3600/np.log(2)
    DipMax=np.min(Dip)

    Scaled=(Dip-DipMax)/(Dip0-DipMax)
    Direct = Dip
    Response = Dip/Dip0

    #plt.plot(LogConc, Direct, "o", lw=2)
    ##plt.plot(LogConc, Response, "o", lw=2, label="dip_rate")

    #####################plt.figure()
    # DipRate###

    # Viability
    ctrl_72h = np.exp(K_Div*72*3600)
    viability = np.exp(np.array(Kprolif)*72*3600)/ctrl_72h
    plt.plot(LogConc, viability, "s", lw=2, label="viability 72h")
    plt.legend(loc=0)



    # CurveFit###     #of dip rate


    def func(x, emax, ec, h):
        eo = 1
        return emax + (eo - emax)/(1+(10**x/ec)**h)

    #this is the dip rate dose response fit
    popt, pcov = curve_fit(func, LogConc, Response, maxfev=50000)
    Emax = popt[0]
    EC = popt[1]
    h = popt[2]
    Eo = 1
    Ei = Emax/Eo + (Eo - Emax)/(Eo*(1+(10**LogConc/EC)**h))           #direct effect
    #Ei = (EC**h)/(EC**h+(10**LogConc)**h)                       #scaled
    #plt.plot(LogConc, Ei, "-", lw = 2)

    print(cell)
    # this is the viability dose response fit
    popt, pcov = curve_fit(func, LogConc, viability, maxfev=500000)
    vEmax = popt[0]
    vEC = popt[1]
    vh = popt[2]
    vEo = 1
    vEi = (vEC ** vh) / (vEC ** vh + (10 ** LogConc) ** vh)
    plt.plot(LogConc, vEi, "-", lw=2)

    plt.xlabel(r'log$_{10}$ conc')
    plt.ylabel("relative effect")

    plt.tight_layout()

    #plt.show()
    #quit()

    #plt.show()############################################################################################################################
    #plt.savefig
    # CurveFit###

    # Debug###
    #print("this is it ",type(Ei))
    #print("this is it ",type(LogConc))
    #print(np.size(Ei))
    #print(np.size(LogConc))
    #print(Ei)
    #print(LogConc)
    #print("this is h ", h)
    #print(EC)
    #print(Emax)
    # Debug###

    # Trapazoid for AUC###
    Eif = np.max(Ei)
    Eis = np.min(Ei)
    Lcf = np.max(LogConc)
    Lcs = np.min(LogConc)

    #print(np.max(Ei),np.min(Ei),np.max(LogConc),np.min(LogConc))

    TrapArea=0
    Count=0
    while Count < (n-1):

        Trap = (Ei[Count]+Ei[Count+1])/2-Eis
        TrapArea=TrapArea+Trap
        Count=Count+1

    print("TrapArea is ",TrapArea)

    # Trapazoid for AUC###

    # Equation Analysis ###
    dmin=10**Lcs
    dmax=10**Lcf

    IC=EC*(1/(1-2*Eis/Eif))**(1/h)

    print("EC is ", EC)
    print("IC is ", IC)
    print("Eif is ", Eif)
    print("Eis is ", Eis)
    print("h is ", h)

    AUC=(1/h)*np.log10(abs((EC**h+dmin**h)/dmin**h*(dmax**h/(EC**h+dmax**h))))
    #AUCA=1/h*np.log10(abs((EC**h+dmin**h)/dmin**h))
    AA=(1/h)*np.log10(abs((EC**h+dmax**h)/(EC**h+dmin**h)))
    #AA=Area-AUC

    print("Area under curve is ",AUC)
    print("Activity Area is ",AA)
    #print (AUCA)
    #plt.show()

    CellAA[cell]=AA
    CellIC[cell]=IC

    # Equation Analysis ###




    # Trapazoid for AUC###
    vEif = np.max(vEi)
    vEis = np.min(vEi)



    dmin=10**Lcs
    dmax=10**Lcf

    vIC=vEC*(1/(1-2*vEis/vEif))**(1/vh)

    print("vEC is ", vEC)
    print("vIC is ", vIC)
    print("vEif is ", vEif)
    print("vEis is ", vEis)
    print("vh is ", vh)
    print("cell is ", cell)

    vAUC=(1/vh)*np.log10(abs((vEC**vh+dmin**vh)/dmin**vh*(dmax**vh/(vEC**vh+dmax**vh))))
    #AUCA=1/h*np.log10(abs((EC**h+dmin**h)/dmin**h))
    vAA=(1/vh)*np.log10(abs((vEC**vh+dmax**vh)/(vEC**vh+dmin**vh)))
    #AA=Area-AUC

    print("vArea under curve is ",vAUC)
    print("vActivity Area is ",vAA)
    #print (AUCA)
    #plt.show()

    vCellAA[cell]=vAUC
    vCellIC[cell]=vIC



    #plt.show()
    #quit()#########################################################################

    # Equation Analysis ###
    cell=cell+1

print("vAA is ", vCellAA)
print("vIC is ", vCellIC)
print("AA is ", CellAA)
print("IC is ", CellIC)

print(type(vCellAA))

vCellAA = np.array(vCellAA)
vCellIC = np.array(vCellIC)
CellAA = np.array(CellAA)
CellIC = np.array(CellIC)

dictionary = {Cell_key[i]:[Cell_div_times[i], vCellAA[i], CellAA[i], vCellIC[i], CellIC[i]] for i in range(len(Cell_key))}
#All_data = dict(zip(Cell_key, zip(Cell_div_times, zip(vCellAA, zip(CellAA, zip(vCellIC, CellIC))))))
#print("All_data is", All_data)
print(dictionary)



# write data to a file
with open('dictionary.json', 'w') as fp:
    json.dump(dictionary, fp)

# read data from a file
with open('dictionary.json', 'r') as fp:
    data = json.load(fp)

print(data)
# Output: {1: [4, 7, 10, 13], 2: [5, 8, 11, 14], 3: [6, 9, 12, 15]}





#All_datas =All_data.sort(key=All_data[2])
#print("All_datas is", All_datas)

plt.figure()
plt.title("Activity Area")
#plt.plot(np.linspace(0, cell, cell),vCellAA,"o", lw = 2, label="vAA")
#plt.plot(np.linspace(0, cell, cell),vCellIC,"o", lw = 2, label="vIC")
#plt.plot(np.linspace(0, cell, cell),CellAA,"o", lw = 2, label="AA")
#plt.plot(np.linspace(0, cell, cell),CellIC,"o", lw = 2, label="IC")
plt.plot(np.linspace(0, 18, 18),vCellAA,"o", lw = 2, label="vAA")
plt.plot(np.linspace(0, 18, 18),CellAA,"o", lw = 2, label="AA")
plt.xlabel('cell line doubling time')
plt.ylabel('value')
plt.legend(loc="best")

plt.savefig("AA_Dict.pdf",format="pdf")


plt.figure()
plt.title("IC_50")
plt.plot(np.linspace(0, 18, 18),np.log10(vCellIC),"o", lw = 2, label="vIC")
plt.plot(np.linspace(0, 18, 18),np.log10(CellIC),"o", lw = 2, label="IC")
plt.xlabel('cell line doubling time')
plt.ylabel('log_10 value')
plt.legend(loc="best")

plt.savefig("IC_Dict.pdf",format="pdf")

plt.figure()
plt.title("Activity Area")
#plt.plot(np.linspace(0, cell, cell),vCellAA,"o", lw = 2, label="vAA")
#plt.plot(np.linspace(0, cell, cell),vCellIC,"o", lw = 2, label="vIC")
#plt.plot(np.linspace(0, cell, cell),CellAA,"o", lw = 2, label="AA")
#plt.plot(np.linspace(0, cell, cell),CellIC,"o", lw = 2, label="IC")
plt.plot(Cell_div_times,vCellAA,"o", lw = 2, label="vAA")
plt.plot(Cell_div_times,CellAA,"o", lw = 2, label="AA")
plt.xlabel('cell line doubling time')
plt.ylabel('value')
plt.legend(loc="best")

plt.savefig("AA_Cell_double.pdf",format="pdf")


plt.figure()
plt.title("IC_50")
plt.plot(Cell_div_times,np.log10(vCellIC),"o", lw = 2, label="vIC")
plt.plot(Cell_div_times,np.log10(CellIC),"o", lw = 2, label="IC")
plt.xlabel('cell line doubling time')
plt.ylabel('log_10 value')
plt.legend(loc="best")

plt.savefig("IC_Cell_double.pdf",format="pdf")



plt.figure()

plt.plot(Cell_div_times,vCellAA-vCellAA[0],"o", lw = 2, label="vAA")
plt.plot(Cell_div_times,np.log10(vCellIC)-np.log10(vCellIC[0]),"o", lw = 2, label="vIC")
plt.plot(Cell_div_times,CellAA-CellAA[0],"o", lw = 2, label="AA")
plt.plot(Cell_div_times,np.log10(CellIC)-np.log10(CellIC[0]),"o", lw = 2, label="IC")
plt.xlabel('cell line')
plt.ylabel('change in value')
plt.legend(loc="best")
plt.savefig("relative_changes.pdf",format="pdf")

plt.show()
