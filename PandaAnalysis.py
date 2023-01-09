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

print(df)
print(df.T)
dfA = df.T

print(dfA)
dfA = dfA.sort_values(by=['vAA'])
sns.heatmap(dfA)
#plt.figure()







plt.show()