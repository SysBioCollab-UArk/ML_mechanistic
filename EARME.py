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


#print (model.observables)


#quit()



#ts = np.linspace(0, 40000, 201)

n = 15

Kprolif = [0]*n
TTD = [0]*n
ProbL = [0]*n
Start = 10**np.linspace(-6, 8, n)
K_Div = 0.017503716680806698/3600
#quit()

TODT= 0
Count = 0
print(Start)


ts = np.linspace(0, 2e5, 1001)

#for init in model.initial_conditions:
 #   print (init)
#quit()


while Count < n:

    sim = ScipyOdeSimulator(model, ts)
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
        print("ERROR: CPARP did not reach 50% of maximum")
        quit()



    ProbL[Count] = Prob

    print("Prob of Life =", format(Prob * 100, ".2f"), "%")
    Kprolif[Count] = K_Div * ProbL[Count] - (1. / TTD[Count]) * (1. - ProbL[Count])




    Count += 1

print("ProbL of Life =", ProbL)
print("Kprolif =", Kprolif)




# DipRate###

for k in Kprolif:
    plt.plot(ts, np.exp(ts*k), label="Kprolif = %g" % k)
plt.xlabel('time')
plt.ylabel('cell count')
plt.legend(loc="best")
plt.yscale("log", base=2)

plt.figure()
Dip = np.array(Kprolif)/np.log(2)
LogConc = np.log10(np.array(Start))

plt.plot(LogConc, Dip/(K_Div/np.log(2)), "o", lw=2)

# DipRate###

# CurveFit###     #of dip rate


def func(x, emax, ec, h):
    e0 = 1
    return emax + (e0 - emax)/(1+(10**x/ec)**h)


popt, pcov = curve_fit(func, LogConc, Dip/(K_Div/np.log(2)), maxfev=50000)
Emax = popt[0]

EC = popt[1]
h = popt[2]

Eo = 1
Ei = Emax + (Eo - Emax)/(1+(10**LogConc/EC)**h)
plt.plot(LogConc, Ei, "-", lw = 2)

plt.xlabel(r'log$_{10}$ conc')
plt.ylabel(r'DIP/DIP$_0$')

plt.show()

# CurveFit###