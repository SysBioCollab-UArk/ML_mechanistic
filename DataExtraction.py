import pandas as pd

# making data frame from csv file
data = pd.read_csv("data/sanger-viability.csv")

# retrieving rows by loc method
data = data.set_index("ARXSPAN_ID")
Cell = data.loc[["ACH-002148"]]
Cell = Cell.set_index("DRUG_NAME")
Drug = Cell.loc[["DOCETAXEL"]]

# display
print(data)
print(Cell)
print(Drug)
print(Drug.to_string())