import pandas as pd

# Load the dataset
df = pd.read_csv("dataset.csv", sep="\t")  # change the separator if necessary

# Remove duplicates on the Smiles column, keeping the first occurrence
df_unique = df.drop_duplicates(subset="Smiles", keep="first")

# Save the new dataset
df_unique.to_csv("dataset_unique.csv", index=False)

print(f"Original rows: {len(df)}")
print(f"Rows after duplicate removal: {len(df_unique)}")
print(f"Duplicates removed: {len(df) - len(df_unique)}")
