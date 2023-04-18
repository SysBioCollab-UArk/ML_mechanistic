import pandas as pd
import seaborn as sns
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

df = pd.read_json('dictionary.json', )
print(df)
df = df.rename(
    index={0: 'Dt', 1: 'vAA', 2: 'AA', 3: 'vIC', 4: 'IC', 5: 'Name'}
)

#df = df[['Name', 'vIC']]



print(df)
print(df.T)
dfA = df.T

dfA=dfA.drop(index=('CAPAN-2'))
dfA=dfA.drop(index=('EFO-21'))
dfA=dfA.drop(index=('RPMI-8226'))



print(dfA)
dfA = dfA[[ 'AA']]
dfA = dfA.sort_values(by=['AA'])
print(dfA)
sns.heatmap(dfA)
#plt.figure()







plt.show()