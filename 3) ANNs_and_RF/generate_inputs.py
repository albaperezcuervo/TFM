# =============================
# Librarys
# =============================
import pandas as pd
import numpy as np
import h5py
# =============================
# CNN Input preparation
# =============================
# Load the OTU taxonomy data
Datdf = pd.read_csv("combined_otu_taxonomy.csv", sep=",")

# Load the metadata file
metadata = pd.read_csv("metadata.csv", sep=",")

# Create a dictionary to map diagnosis categories to binary values
diagnosis_map = {"Control": 0, "Alzheimer": 1}

# Map diagnosis values to binary and add as a new column in the metadata
metadata["Diagnosis_Binary"] = metadata["diagnosis"].map(diagnosis_map)

# Extract columns starting with "SRR" from the primary DataFrame `df`
srr_columns = [col for col in df.columns if col.startswith("SRR")]

# Build a DataFrame mapping each SRR column to its corresponding binary diagnosis
ordered_diagnosis_df = (
    pd.DataFrame({"SRR_Column": srr_columns})
    .merge(metadata[["SampleID", "Diagnosis_Binary"]], 
           left_on="SRR_Column", 
           right_on="SampleID", 
           how="left")
    .drop(columns=["SampleID"])  # Drop redundant SampleID column
)

# Display the intermediate mapping DataFrame
print(ordered_diagnosis_df)

# Convert the mapping DataFrame to a 2D NumPy array
srr_diagnosis_array = ordered_diagnosis_df[["SRR_Column", "Diagnosis_Binary"]].to_numpy()

# Display the resulting array
print(srr_diagnosis_array)

# Save the array to a .npy file
np.save('metadata_array.npy', srr_diagnosis_array)

# Re-extract SRR columns from the main DataFrame
srr_columns = [col for col in df.columns if col.startswith('SRR')]

# Initialize a list to store individual arrays for each SRR sample
rd_arrays = []

# Iterate over each SRR column to process and extract relevant features
for srr in srr_columns:
    # Create temporary DataFrame with new feature columns initialized to zero
    new_columns = pd.DataFrame({
        f'k_{srr}': [0] * len(df),
        f'p_{srr}': [0] * len(df),
        f'c_{srr}': [0] * len(df),
        f'o_{srr}': [0] * len(df),
        f'f_{srr}': [0] * len(df),
        f'g_{srr}': df[srr]
    })
    
    # Append the new feature columns to the main DataFrame
    df = pd.concat([df, new_columns], axis=1)
    
    # Assign 1 when a change in "Kingdom" is detected
    df[f'k_{srr}'] = (df['Kingdom'] != df['Kingdom'].shift()).astype(int)

    # Compute accumulated abundance for each taxonomic level
    for level, col_prefix in zip(['Family', 'Order', 'Class', 'Phylum'], ['f_', 'o_', 'c_', 'p_']):
        df[f'{col_prefix}{srr}'] = df.groupby(level)[f'g_{srr}'].transform('sum')
        df[f'{col_prefix}{srr}'] = df[level].ne(df[level].shift()).astype(int) * df[f'{col_prefix}{srr}']

    # Select the columns of interest and convert them to a NumPy array
    filtered_array = df[[f'k_{srr}', f'p_{srr}', f'c_{srr}', f'o_{srr}', f'f_{srr}', f'g_{srr}']].to_numpy()

    # Append the processed array to the list
    rd_arrays.append(filtered_array)

# Convert the list of arrays into a single NumPy array
main_array = np.array(rd_arrays)

# Display the shape of the resulting array
print(f"Global array dimensions: {main_array.shape}")

# Add a final dimension to match the required shape (samples, rows, features, 1)
main_array = main_array[..., np.newaxis]

# Confirm the new shape
print("Shape of main_array:", main_array.shape)  # Expected: (60, 248, 6, 1)

# Save the final array to an HDF5 file
with h5py.File('merged_abundance_acumulative.h5', 'w') as hf:
    hf.create_dataset('dataset_name', data=main_array)


# =============================
# MLPNNs and RF input preparation
# =============================

# Extract SRR columns from the previously loaded DataFrame
srr_columns = df.filter(like='SRR').columns

# Convert each SRR column into a separate NumPy array
numpy_arrays = [df[srr].to_numpy() for srr in srr_columns]

# Combine the list of arrays into a single NumPy array
numpy_array_final = np.array(numpy_arrays)

# Display the final NumPy array for verification
print(numpy_array_final)

# Save the NumPy array to an HDF5 file without flattening
with h5py.File('Merge_Abundance_genus.h5', 'w') as hf:
    hf.create_dataset('dataset_name', data=numpy_array_final)
