import os
from Bio.PDB import PDBParser, Select, PDBIO
import Bio.PDB as pdb
complex_pdb = r"C:\Users\123\Desktop\jobwork\subject\IL-8\结构\colab\SMAFHSH.pdb"
monomers_dir = r"C:\Users\123\Desktop\jobwork\subject\IL-8\结构\colab\SMAFHSH"


# 打开复合体PDB文件  
parser = pdb.PDBParser()
structure = parser.get_structure('complex', complex_pdb)

# 遍历所有模型
for model in structure:
    # 遍历所有链 
    for chain in model:
        # 输出单体PDB 
        output_filename = f'{chain.id}.pdb' # 定义单体PDB文件名
        io = pdb.PDBIO()
        io.set_structure(chain)
        io.save(os.path.join(monomers_dir, output_filename))
    chains = ['A','B']
    IL8_chain = pdb.Model.Model(0)
    for chain in chains:
        if chain in chains:
            IL8_chain.add(chain)  

    # 输出整个复合体PDB
    output_filename = 'IL8.pdb'
    io = pdb.PDBIO()
    io.set_structure(IL8_chain) 
    io.save(os.path.join(monomers_dir, output_filename))