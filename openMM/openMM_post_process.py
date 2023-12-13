import pandas as pd
import ast
import csv

# Assuming you have already read the CSV file and created the original DataFrame df1

with open (csv_path) as f:
    reader = csv.reader(f)
    for row in reader:
        data = list(reader)
        df = pd.DataFrame(data)
        df1 = df.T

# Function to parse the custom string
def parse_energy_string(s):
    # Extracting information using string manipulation
    file_name_start = s.find("file_name='") + len("file_name='")
    file_name_end = s.find("'", file_name_start)
    file_name = (s[file_name_start:file_name_end]).split('_')[0]

    raw_start = s.find("raw=EnergyUnit(") + len("raw=EnergyUnit(")
    raw_end = s.find(")", raw_start)
    raw_data = s[raw_start:raw_end].split(', ')

    raw_complex_energy = float(raw_data[0].split('=')[1])
    raw_protein_energy = float(raw_data[1].split('=')[1])
    raw_peptide_energy = float(raw_data[2].split('=')[1])
    raw_diff_energy = float(raw_data[3].split('=')[1])

    minimized_start = s.find("minimized=EnergyUnit(") + len("minimized=EnergyUnit(")
    minimized_end = s.find(")", minimized_start)
    minimized_data = s[minimized_start:minimized_end].split(', ')

    minimized_complex_energy = float(minimized_data[0].split('=')[1])
    minimized_protein_energy = float(minimized_data[1].split('=')[1])
    minimized_peptide_energy = float(minimized_data[2].split('=')[1])
    minimized_diff_energy = float(minimized_data[3].split('=')[1])

    return pd.Series([file_name, raw_complex_energy, raw_protein_energy, raw_peptide_energy, raw_diff_energy,
                      minimized_complex_energy, minimized_protein_energy, minimized_peptide_energy, minimized_diff_energy])

# Apply the function to the first column of the DataFrame
new_df = df1[0].apply(parse_energy_string)

# Concatenate the new DataFrame with the original DataFrame
result_df = pd.concat([df1, new_df], axis=1)

# Set appropriate column names
result_df.columns = ['Original_Column'] + ['file_name', 'raw_complex_energy', 'raw_protein_energy', 'raw_peptide_energy', 'raw_diff_energy',
                                           'minimized_complex_energy', 'minimized_protein_energy', 'minimized_peptide_energy', 'minimized_diff_energy']

# Drop the original column if you don't need it anymore
result_df = result_df.drop(columns=['Original_Column'])

result_df.to_csv('output_file.csv', index=False)
# Display the resulting DataFrame
print(result_df)

