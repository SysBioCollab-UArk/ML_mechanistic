import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import numpy as np
import seaborn as sns

# cell_line Key
# 0      A549       ACH-000681
# 1      AN3-CA     ACH-000940
# 2      CAPAN-2    ACH-000107
# 3      Daudi      ACH-000786
# 4      EFO-21     ACH-000308
# 5      HDLM-2     ACH-000267
# 6      HL-60      ACH-000002
# 7      HeLa       ACH-001086
# 8      K-562      ACH-000551
# 9      MCF7       ACH-000019
# 10     MOLT-4     ACH-001127
# 11     PC-3       ACH-000090
# 12     REH        ACH-000960
# 13     RPMI-8226  ACH-000817
# 14     RT4        ACH-000242
# 15     SiHa       ACH-000556
# 16     THP-1      ACH-000146
# 17     WM-115     ACH-000304

# enter for CellKeyNum, the cell_line number
CellKeyNum = 0
# enter for s, the smoothing factor (running average frame)
s = 0
# enter for so, 1 for viability to be driven to zero
so = 1
# enter for Dname, the drug name

# Dname = "DOCETAXEL"
Dname = "RTRAIL"

#Dname = "NAVITOCLAX"
#Dname = "OBATOCLAX"
#Dname = "TW-37"

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

print("new_dict", new_dict)

# print(test)

values_list = list(new_dict.values())
keys_list = list(new_dict.keys())
NameR = {k: v for k, v in zip(values_list, keys_list)}


print(values_list)

# making data frame from csv file
data = pd.read_csv("data/sanger-dose-response.csv")
print(data.columns)
data = data[['DRUG_NAME', 'DRUG_ID', 'COSMIC_ID', 'DATASET', 'ARXSPAN_ID',
             'IC50_PUBLISHED','AUC_PUBLISHED','ec50','upper_limit','lower_limit','mse'




]]

#


#quit()





filtered_df = data.loc[data['ARXSPAN_ID'].isin(values_list[CellKeyNum:CellKeyNum + 18])]
# retrieving rows by loc method
# data = data.set_index("ARXSPAN_ID")
Cell = filtered_df  # data  # filtered_df
# Cell = data.loc[["ACH-002148"]]
Cell = Cell.set_index("DRUG_NAME")
Drug = Cell.loc[[Dname]]
Drug = Drug.set_index("DATASET")
# Drug = Drug.loc[["GDSC1"]]
#Drug.sort_values(by='dose', ascending=True, inplace=True)
# display
print(data)
print(Cell)
print("start",Drug,"end")
#print(Drug.to_string())


Drug.to_csv('my_data.csv', index=False)

Drug = Drug[[ 'ARXSPAN_ID',
             'AUC_PUBLISHED','mse']]








print(NameR)

#quit()

Drug['ARXSPAN_ID'] = Drug['ARXSPAN_ID'].map(NameR).fillna(Drug['ARXSPAN_ID'])

print("start",Drug,"end")


Drug.to_csv('my_data.csv', index=False)

Drug = Drug.sort_values(by=['AUC_PUBLISHED'])
sns.heatmap(Drug.pivot_table(values='AUC_PUBLISHED', index='ARXSPAN_ID'))



plt.show()