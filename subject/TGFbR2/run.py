import os 
import csv
import sys 
import pandas as pd
sys.path.append(os.path.abspath(''))
from openMM.caculate_potential_energy import batch_calc_and_rank_by_raw_diff_energy

def run(data_path):
    cache_csv = os.path.join(data_path, "cache.csv")
    pdb_paths = []
    for root, dirs, files in os.walk(data_path):
        for file in files:
            if file.endswith(".pdb"):
                pdb_path = os.path.join(root, file)
                pdb_paths.append(pdb_path)
    sorted_index, sorted_energies = batch_calc_and_rank_by_raw_diff_energy(pdb_paths)
                #save sorted_index, sorted_energies to csv        
    with open(cache_csv, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(sorted_index)
        writer.writerow(sorted_energies)

    return cache_csv

def read_energy_csv(csv_path):
    with open (csv_path) as f:
        reader = csv.reader(f)
        for row in reader:
            data = list(reader)
            df = pd.DataFrame(data)
            df1 = df.T
    return df1

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

def result_csv(data_path):
    #cache_csv = run(data_path)
    df = read_energy_csv(cache_csv)

    new_df = df[0].apply(parse_energy_string)
    result_df = pd.concat([df, new_df], axis=1)
    result_df.columns = ['Original_Column'] + ['file_name', 'raw_complex_energy', 'raw_protein_energy', 'raw_peptide_energy', 'raw_diff_energy',
                                            'minimized_complex_energy', 'minimized_protein_energy', 'minimized_peptide_energy', 'minimized_diff_energy']

    # Drop the original column if you don't need it anymore
    result_df = result_df.drop(columns=['Original_Column'])
    result_csv_path = os.path.join(data_path, 'result.csv')
    #save sorted_index, sorted_energies to csv
    result_df.to_csv(result_csv_path, index=False)



if __name__ == "__main__":
    data_path = '/mnt/nas1/lanwei-125/TGFbR2/CF_relax/extracted_files/'
    cache_csv = '/mnt/nas1/lanwei-125/TGFbR2/CF_relax/extracted_files/cache.csv'
    result_csv(data_path)