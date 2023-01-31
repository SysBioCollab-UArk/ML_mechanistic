import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# Read in the "doubling_times.csv" file
df = pd.read_csv("data/doubling_times.csv")

# Filter the "cell_line" and "depmap_id" columns into a new dataframe
cell_df = df[['cell_line', 'depmap_id']].drop_duplicates().set_index("cell_line")

# Convert the dataframe into a dictionary with cell_line as key and depmap_id as value
cell_dict = cell_df.to_dict("split")['data']
cell_dict = {k: v[0] for k, v in cell_dict.items()}

# Read in the "sanger-viability.csv" file
data = pd.read_csv("data/sanger-viability.csv")

# Filter the columns to include "DRUG_NAME","dose","viability","ARXSPAN_ID"
data = data[['DRUG_NAME','dose','viability','ARXSPAN_ID']]

# Filter rows where ARXSPAN_ID is in the cell_dict
data = data[data['ARXSPAN_ID'].isin(cell_dict.values())]

# Set the index to "DRUG_NAME"
data = data.set_index("DRUG_NAME")

# Filter only rows with "DOCETAXEL"
docetaxel_df = data.loc[["DOCETAXEL"]]

# Plot line and scatter plots of the viability against dose
docetaxel_df.plot(x='dose', y='viability', kind='line')
docetaxel_df.plot.scatter(x='dose', y='viability')

# Plot bar plot of mean viability grouped by DRUG_NAME
docetaxel_df.groupby('DRUG_NAME')['viability'].mean().plot(kind='bar')

# Plot heatmap of pivot table of viability indexed by DRUG_NAME and columns dose
sns.heatmap(docetaxel_df.pivot_table(values='viability', index='DRUG_NAME', columns='dose'))

# Show all plots
plt.show()

