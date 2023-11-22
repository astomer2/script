
#将脚本拖进pymol的命令行里，然后输入run excute_ADCP_PDB.py
#或者在pymol的命令行里输入run C:\Users\123\Desktop\jobwork\subject\script\excute_ADCP_PDB.py
#最后输入excute_ADCP_PDB path/name_list,path/protein, path
from pymol import cmd
import os
import glob

def excute_ADPC_PDB(protein ,path):
    pdb_file_name = [os.path.basename(x) for x in glob.glob(path + '*.pdb')]
    #name_list = os.path.join(path, "name.txt")
    #with open(name_list) as f:
    for line in pdb_file_name:
        line = line.strip()
        peptide = f'{path}/{line}.pdb'
        cmd.load(protein)
        cmd.load(peptide)
        cmd.alter(line,'chain="C"')
        cmd.save(peptide)
        cmd.delete(all)
        print(peptide)

cmd.extend('excute_ADPC_PDB', excute_ADPC_PDB)


