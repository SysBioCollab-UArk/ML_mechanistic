import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import numpy as np
import seaborn as sns

#cell_line Key
#0      A549       ACH-000681
#1      AN3-CA     ACH-000940
#2      CAPAN-2    ACH-000107
#3      Daudi      ACH-000786
#4      EFO-21     ACH-000308
#5      HDLM-2     ACH-000267
#6      HL-60      ACH-000002
#7      HeLa       ACH-001086
#8      K-562      ACH-000551
#9      MCF7       ACH-000019
#10     MOLT-4     ACH-001127
#11     PC-3       ACH-000090
#12     REH        ACH-000960
#13     RPMI-8226  ACH-000817
#14     RT4        ACH-000242
#15     SiHa       ACH-000556
#16     THP-1      ACH-000146
#17     WM-115     ACH-000304

# enter for x the cell_line number
x=0



name = pd.read_csv("data/doubling_times.csv")

print(name.to_string())
Name = name[['cell_line', 'depmap_id']]
Name = Name.set_index("cell_line")
Name = Name.drop_duplicates()
print(Name)


test = Name.to_dict("split")
data = test
new_dict = {}
for i in range(len(data['index'])):
    new_dict[data['index'][i]] = data['data'][i][0]

print("new_dict",new_dict)


#print(test)

values_list = list(new_dict.values())
keys_list = list(new_dict.keys())
print(values_list)


# making data frame from csv file
data = pd.read_csv("data/sanger-viability.csv")
data = data[['DRUG_NAME','DRUG_ID','COSMIC_ID','dose','viability','DATASET','ARXSPAN_ID']]


filtered_df = data.loc[data['ARXSPAN_ID'].isin(values_list [x:x+1])]




# retrieving rows by loc method
#data = data.set_index("ARXSPAN_ID")


Dname = "NAVITOCLAX"
#Dname = "DOCETAXEL"



Cell = filtered_df
#Cell = data.loc[["ACH-002148"]]
Cell = Cell.set_index("DRUG_NAME")
Drug = Cell.loc[[Dname]]
Drug = Drug.set_index("DATASET")
#Drug = Drug.loc[["GDSC1"]]
Drug.sort_values(by='dose', ascending=True, inplace=True)








# display
print(data)
print(Cell)
print(Drug)
print(Drug.to_string())


# from data frame to numpy array
viability = Drug['viability'].to_numpy()
Conc = Drug['dose'].to_numpy()
LogConc = np.log10(Conc)

# viability plot
plt.plot(LogConc, viability, "s", lw=2, label="viability")
plt.legend(loc=0)

# Curve fit viability

def func(x, emax, ec, h):
    eo = 1
    return emax + (eo - emax) / (1 + (10 ** x / ec) ** h)



popt, pcov = curve_fit(func, LogConc, viability, maxfev=500000)
vEmax = popt[0]
vEC = popt[1]
vh = popt[2]
vEo = 1
vEi = (vEC ** vh) / (vEC ** vh + (10 ** LogConc) ** vh)
plt.plot(LogConc, vEi, "-", lw=2)
plt.xlabel(r'log$_{10}$ conc')
plt.ylabel("relative effect")
plt.title(Dname + ", " + keys_list[x])
plt.tight_layout()


# variables
vEif = np.max(vEi)
vEis = np.min(vEi)
Lcf = np.max(LogConc)
Lcs = np.min(LogConc)

# area under curve
dmin=10**Lcs
dmax=10**Lcf

vIC=vEC*(1/(1-2*vEis/vEif))**(1/vh)

vAUC=(1/vh)*np.log10(abs((vEC**vh+dmin**vh)/dmin**vh*(dmax**vh/(vEC**vh+dmax**vh))))
vAA=(1/vh)*np.log10(abs((vEC**vh+dmax**vh)/(vEC**vh+dmin**vh)))

print("vArea under curve is ",vAUC)
print("vActivity Area is ",vAA)
print("vEC is ", vEC)
print("vIC is ", vIC)



# various plots

#Drug.plot(x='dose', y='viability', kind='line')
#plt.figure()
Drug.plot.scatter(x='dose', y='viability')

#Drug.groupby(['dose','DATASET'])['viability'].mean().plot(kind='bar')
#plt.figure()
#sns.heatmap(Drug.pivot_table(values='viability', index='DRUG_NAME', columns='dose'))
plt.show()