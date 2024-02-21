from pathlib import Path
import os 
import pandas as pd
import matplotlib.pyplot as plt

def get_residue_num(pdb_path):
    """
    Counts the number of unique residues in a PDB file.

    Args:
        pdb_path (str): The path to the PDB file.

    Returns:
        int: The number of unique residues in the PDB file.
    """
    with open(pdb_path) as f:
        lines = f.readlines()

    residue_count = 0
    prev_res_id = 0
    peptide_id = 0 
    for line in lines:
        if line.startswith('ATOM'):
            peptide_id = line[22:26].strip()
            if peptide_id != prev_res_id:
                residue_count += 1
            prev_res_id = peptide_id
    return residue_count

root_dir = Path('/mnt/nas1/lanwei-125/MC5R/MD/pos/ZMK2-110/')

file_path = root_dir / 'decomp-energy.csv'
data = pd.read_csv(file_path, skiprows=8, delimiter=',',header=None, engine='python')

a = data[0].str.split(expand=True, n=1)
data['Residue ID'] = a[0]
data['ids'] = a[1]

# 将第 1 列拆分为三个列
residue_info = data[1].str.split(expand=True, n=2)
data['chain ID'] = residue_info[0]
data['Amino Acid'] = residue_info[1]
data['Residue_ID'] = residue_info[2]

# 删除第 0 列和第 1 列
#data = data.drop([0, 1], axis=1)
data["energy"] = data[17]

A= get_residue_num(root_dir/'protein.pdb')
data_subset = data.iloc[:A]
residue_ids = data_subset['Residue_ID']
three_code_id = data_subset[0]
energy_contributions = data_subset['energy']



plt.plot(residue_ids, energy_contributions, label='a-MSH')
plt.title('Energy Contributions of Amino Acid Residues')
plt.xlabel('Residue ID')
plt.xlabel
plt.ylabel('Energy Contribution (Kcal/mol)')
plt.legend()
plt.xticks(ticks= range(-1,A,20))
for x, y, three_code_id in zip(residue_ids, energy_contributions, three_code_id):
    if abs(y)>0.1:
        print(three_code_id)
        plt.annotate(f"{three_code_id}", (x, y), xytext=(5, 5), textcoords='offset points', ha='left'
                     ,arrowprops=dict(arrowstyle='->') 
                     )

plt.tight_layout() 
plt.gcf().set_size_inches(20, 8)


#plt.savefig(root_dir / 'MC5R_energy_contributions.png',dpi=600)
plt.show()

