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

import read_data
from read_data import get_k_div

Cell_div_times=list(read_data.get_k_div("data/doubling_times.csv").values())
Cell_key=list(read_data.get_k_div("data/doubling_times.csv").keys())

print(Cell_key)
print(Cell_div_times)
#quit()



Cell_div_times=np.array([0.02469112365984261,0.013524823035316005,0.009973340727481227, 0.01933464938800405, 0.015322684707684294, 0.009560650766344074,0.024755256448569473,0.017503716680806698,0.0204669443078724,0.016681350442654508,0.021061901566695386,0.02161132344294155,0.016675553662228032 ,0.013119504994825462,0.011375500775053806,0.011108127893588867,0.015576341136178546,0.018216745875425634])

CellAA = [0]*len(Cell_div_times)
CellIC = [0]*len(Cell_div_times)
vCellAA = [0]*len(Cell_div_times)
vCellIC = [0]*len(Cell_div_times)
ts = np.linspace(0, 72*60*60, 1001)
cell = 0
sim = ScipyOdeSimulator(model, ts)
while cell < len(Cell_div_times):

    #print (model.observables)


    n = 30

    Kprolif = [0]*n
    TTD = [0]*n
    ProbL = [0]*n
    Start = 10**np.linspace(-15, 8, n)
    K_Div = Cell_div_times[cell]/3600
    #quit()

    TODT= 0
    Count = 0
    print(Start)


    while Count < n:


        traj = sim.run(initials={model.monomers["L"](b=None): Start[Count]})

        ######################################################################################################
        #print(traj.dataframe["CPARP_total"])

        #plt.plot(traj.dataframe["CPARP_total"])

        #plt.show()
       # quit()
        P = 0
        Q = 0
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
        TOD = D * 200
        print("TOD", TOD)

        Ya = D
        Yb = D + 1
        Xa = TODF
        Xb = TODT
        X = model.parameters["PARP_0"].value/2
        Y = 1
        Y = Ya + (X - Xa) * ((Yb - Ya) / (Xb - Xa))
        print("Y =", Y)
        TOD = Y * 200

        ########################################################################################################

        TTD[Count] = TOD
        Death = TOD     # Time of death
        domain = ts        # timeScale
        x = domain

        ##############################################################################################################

        Prob = 1-np.exp(-K_Div*Death)

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
    print("Kprolif =", Kprolif)

    # DipRate###
    plt.subplot(1, 2, 1)
    plt.title("linear_tc")
    plt.plot(ts/3600, np.exp(ts * K_Div), label="Kprolif = %g" % K_Div)
    for k in Kprolif:
        plt.plot(ts/3600, np.exp(ts * k), label="Kprolif = %g" % k)
    plt.xlabel('time (hr)')
    plt.ylabel('cell count')
    #plt.legend(loc="best")

    plt.subplot(1, 2, 2)
    plt.title("log_tc")
    plt.plot(ts/3600,np.log2(np.exp(ts * K_Div)), label="Kprolif = %g" % K_Div)
    for k in Kprolif:
        plt.plot(ts/3600, np.log2(np.exp(ts*k)), label="Kprolif = %g" % k)
    plt.xlabel('time (hr)')
    plt.ylabel('cell count')
    #plt.legend(loc="best")

    #plt.yscale("log", base=2)

    plt.figure()
    Dip = np.array(Kprolif)/np.log(2)
    LogConc = np.log10(np.array(Start))

    Direct = Dip/(K_Div/np.log(2))
    Scaled = (Direct-np.min(Direct))/(np.max(Direct)-np.min(Direct))

    #plt.plot(LogConc, Direct, "o", lw=2)
    plt.plot(LogConc, Scaled, "o", lw=2, label="dip_rate")


    # DipRate###

    # Viability
    ctrl_72h = np.exp(K_Div*72*3600)
    viability = np.exp(np.array(Kprolif)*72*3600)/ctrl_72h
    plt.plot(LogConc, viability, "s", lw=2, label="viability 72h")
    plt.legend(loc=0)




    # CurveFit###     #of dip rate


    def func(x, emax, ec, h):
        e0 = 1
        return emax + (e0 - emax)/(1+(10**x/ec)**h)

    #this is the dip rate dose response fit
    popt, pcov = curve_fit(func, LogConc, Dip/(K_Div/np.log(2)), maxfev=50000)
    Emax = popt[0]
    EC = popt[1]
    h = popt[2]
    Eo = 1
    #Ei = Emax + (Eo - Emax)/(1+(10**LogConc/EC)**h)            #direct effect
    Ei = (EC**h)/(EC**h+(10**LogConc)**h)                       #scaled
    plt.plot(LogConc, Ei, "-", lw = 2)


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

    vAUC=(1/vh)*np.log10(abs((vEC**vh+dmin**vh)/dmin**vh*(dmax**vh/(vEC**vh+dmax**vh))))
    #AUCA=1/h*np.log10(abs((EC**h+dmin**h)/dmin**h))
    vAA=(1/vh)*np.log10(abs((vEC**vh+dmax**vh)/(vEC**vh+dmin**vh)))
    #AA=Area-AUC

    print("vArea under curve is ",vAUC)
    print("vActivity Area is ",vAA)
    #print (AUCA)
    #plt.show()

    vCellAA[cell]=vAA
    vCellIC[cell]=vIC

    # Equation Analysis ###
    cell=cell+1
print("vAA is ", vCellAA)
print("vIC is ", vCellIC)
print("AA is ", CellAA)
print("IC is ", CellIC)


vCellAA = [x - vCellAA[0] for x in vCellAA]
vCellIC = [x - vCellIC[0] for x in vCellIC]
CellAA = [x - CellAA[0] for x in CellAA]
CellIC = [x - CellIC[0] for x in CellIC]

plt.figure()
plt.title("Standardized Metrics")
#plt.plot(np.linspace(0, cell, cell),vCellAA,"o", lw = 2, label="vAA")
#plt.plot(np.linspace(0, cell, cell),vCellIC,"o", lw = 2, label="vIC")
#plt.plot(np.linspace(0, cell, cell),CellAA,"o", lw = 2, label="AA")
#plt.plot(np.linspace(0, cell, cell),CellIC,"o", lw = 2, label="IC")
plt.plot(Cell_div_times,vCellAA,"o", lw = 2, label="vAA")
plt.plot(Cell_div_times,vCellIC,"o", lw = 2, label="vIC")
plt.plot(Cell_div_times,CellAA,"o", lw = 2, label="AA")
plt.plot(Cell_div_times,CellIC,"o", lw = 2, label="IC")
plt.xlabel('cell line')
plt.ylabel('change in value')
plt.legend(loc="best")


plt.show()
