from scipy.stats import norm
from scipy.special import erf

from earm.lopez_embedded import model
import numpy as np
import matplotlib.pyplot as plt
from pysb.simulator import ScipyOdeSimulator
from pysb import *
import math
from scipy.optimize import curve_fit

ts = np.linspace(0, 40000, 201)

n = 15

Kprolif = [0]*n
TTD = [0]*n
ProbL = [0]*n
Start = 10**np.linspace(1, 7, n)
Count = 0
print(Start)
while Count < n:

    sim = ScipyOdeSimulator(model, ts)
    traj = sim.run(initials={model.monomers["L"](bf=None): Start[Count]})

    ######################################################################################################
    # print(traj.dataframe["cPARP"])
    P = 0
    Q = 0
    for T in traj.dataframe["cPARP"]:
        # print(T)
        # print(P)
        P = P + 1
        if T < 500000:
            TODF = T
            D = P - 1
        if T > 500000 and Q == 0:
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
    X = 500000
    Y = 1
    Y = Ya + (X - Xa) * ((Yb - Ya) / (Xb - Xa))
    print("Y =", Y)
    TOD = Y * 200

    ########################################################################################################

    TTD[Count] = TOD
    Death = TOD     # Time of death
    Life = 3000   # Avg time cell reproduction
    LDis = 1000      # standard distrabution of Life

    domain = ts        # timeScale
    x = domain

    # ###########Z-Score Calculation#######################

    Zscore = ((Death-Life)/LDis)

    # print("Zscore =", Zscore)

    ####

    w = np.sign(abs(Death-Life)/((2**.5)*LDis))
    q = erf((Death-Life)/((2**.5)*LDis))

    Prob = .5*(1+w * q)

    print("Prob of Life =", format(Prob*100, ".2f"), "%")

    ProbL[Count] = Prob
    Kprolif[Count] = (1 / Life) * ProbL[Count] - (1. / TTD[Count]) * (1. - ProbL[Count])
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

plt.plot(LogConc, Dip/Dip[0], "o-", lw=2)

# DipRate###
# plt.show()
# quit()
# CurveFit###     #of dip rate


def func(x, Emax, EC, h):
    return Emax + (1 - Emax)/(1+(10**x/EC)**h)


popt, pcov = curve_fit(func, LogConc, Dip/Dip[0], maxfev=50000)
Emax = popt[0]

EC = popt[1]
h = popt[2]

Eo = 1
Ei = Emax + (Eo - Emax)/(1+(10**LogConc/EC)**h)
plt.plot(LogConc, Ei, "s")

plt.xlabel(r'log$_{10}$ conc')
plt.ylabel(r'DIP/DIP$_0$')

plt.show()

# CurveFit###

# quit()

# BellCurve###
# plt.figure()
#
# NORM = norm.pdf(domain,Life,LDis)   # normal distrabution with mean=first number standared deviation=second
#
# plt.plot(domain, (NORM * 1000000000))
#
# plt.plot(traj.dataframe["cPARP"], label='cPARP')
# plt.title('standard normal')
# plt.xlabel('value')
# plt.ylabel('Density')
# plt.axvline(Death, color='g')
# plt.axhline(0, color="black")
# plt.fill_between(domain, NORM * 1000000000, where = [(x > 0) and (x <= Death) for x in x] ,color="r")
#
#
# plt.show()
# BellCurve###
