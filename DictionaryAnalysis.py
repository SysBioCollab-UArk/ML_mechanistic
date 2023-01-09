
import numpy as np
import matplotlib.pyplot as plt
from pysb.simulator import ScipyOdeSimulator
from pysb import *
import math
from scipy.optimize import curve_fit
import os
import read_data
import json




with open('dictionary.json', 'r') as fp:
    data = json.load(fp)

print(data)


dictionary = data

AA =[x[2] for x in dictionary.values()]
IC= [x[4] for x in dictionary.values()]
vAA= [x[1] for x in dictionary.values()]
vIC= [x[3] for x in dictionary.values()]




dictionary = {k: v for k, v in zip(AA, vAA)}
sorted_dictionary = {k: v for k, v in sorted(dictionary.items(), key=lambda item: item[0])}

dictionaryA = {k: v for k, v in zip(data.keys(), AA)}
sorted_dictionaryA = {k: v for k, v in sorted(dictionaryA.items(), key=lambda item: item[1])}
print(sorted_dictionaryA)

print(sorted_dictionary)
print(sorted_dictionaryA.keys())




#x_values = list(dictionary.keys())
space = len(list(sorted_dictionaryA.keys()))
print(space)
#x_values = np.linspace(0,space+.5,space)
x_values = np.arange(0,space,1)
y_valuesA = list(sorted_dictionary.keys())
y_valuesB = list(sorted_dictionary.values())
plt.plot(x_values, y_valuesA,"o", lw = 2, label="AA")
plt.plot(x_values, y_valuesB,"o", lw = 2, label="vAA")
plt.legend(loc="best")
plt.xticks(x_values, list(sorted_dictionaryA.keys()))

plt.show()



#dictionary = {Cell_key[i]:[Cell_div_times[i], vCellAA[i], CellAA[i], vCellIC[i], CellIC[i]] for i in range(len(Cell_key))}


#plt.plot(np.linspace(0, cell, cell),vCellAA,"o", lw = 2, label="vAA")










#
#plt.figure()
#plt.title("Activity Area")
##plt.plot(np.linspace(0, cell, cell),vCellAA,"o", lw = 2, label="vAA")
##plt.plot(np.linspace(0, cell, cell),vCellIC,"o", lw = 2, label="vIC")
##plt.plot(np.linspace(0, cell, cell),CellAA,"o", lw = 2, label="AA")
##plt.plot(np.linspace(0, cell, cell),CellIC,"o", lw = 2, label="IC")
#plt.plot(np.linspace(0, 18, 18),vCellAA,"o", lw = 2, label="vAA")
#plt.plot(np.linspace(0, 18, 18),CellAA,"o", lw = 2, label="AA")
#plt.xlabel('cell line doubling time')
#plt.ylabel('value')
#plt.legend(loc="best")
#
#plt.savefig("AA_Dict.pdf",format="pdf")
#
#
#plt.figure()
#plt.title("IC_50")
#plt.plot(np.linspace(0, 18, 18),np.log10(vCellIC),"o", lw = 2, label="vIC")
#plt.plot(np.linspace(0, 18, 18),np.log10(CellIC),"o", lw = 2, label="IC")
#plt.xlabel('cell line doubling time')
#plt.ylabel('log_10 value')
#plt.legend(loc="best")
#
#plt.savefig("IC_Dict.pdf",format="pdf")
#
#plt.figure()
#plt.title("Activity Area")
##plt.plot(np.linspace(0, cell, cell),vCellAA,"o", lw = 2, label="vAA")
##plt.plot(np.linspace(0, cell, cell),vCellIC,"o", lw = 2, label="vIC")
##plt.plot(np.linspace(0, cell, cell),CellAA,"o", lw = 2, label="AA")
##plt.plot(np.linspace(0, cell, cell),CellIC,"o", lw = 2, label="IC")
#plt.plot(Cell_div_times,vCellAA,"o", lw = 2, label="vAA")
#plt.plot(Cell_div_times,CellAA,"o", lw = 2, label="AA")
#plt.xlabel('cell line doubling time')
#plt.ylabel('value')
#plt.legend(loc="best")
#
#plt.savefig("AA_Cell_double.pdf",format="pdf")
#
#
#plt.figure()
#plt.title("IC_50")
#plt.plot(Cell_div_times,np.log10(vCellIC),"o", lw = 2, label="vIC")
#plt.plot(Cell_div_times,np.log10(CellIC),"o", lw = 2, label="IC")
#plt.xlabel('cell line doubling time')
#plt.ylabel('log_10 value')
#plt.legend(loc="best")
#
#plt.savefig("IC_Cell_double.pdf",format="pdf")
#
#
#
#plt.figure()
#
#plt.plot(Cell_div_times,vCellAA-vCellAA[0],"o", lw = 2, label="vAA")
#plt.plot(Cell_div_times,np.log10(vCellIC)-np.log10(vCellIC[0]),"o", lw = 2, label="vIC")
#plt.plot(Cell_div_times,CellAA-CellAA[0],"o", lw = 2, label="AA")
#plt.plot(Cell_div_times,np.log10(CellIC)-np.log10(CellIC[0]),"o", lw = 2, label="IC")
#plt.xlabel('cell line')
#plt.ylabel('change in value')
#plt.legend(loc="best")
#plt.savefig("relative_changes.pdf",format="pdf")
#
#plt.show()
#