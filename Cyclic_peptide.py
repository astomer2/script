import os
import re
from pathlib import Path
import matplotlib.pyplot as plt

residues_dict = {
  'GLY':'G',
  'ALA':'A',
  'SER':'S',
  'THR':'T',
  'CYS':'C',
  'VAL':'V',
  'LEU':'L',
  'ILE':'I',
  'MET':'M',
  'PRO':'P',
  'PHE':'F',
  'TYR':'Y',
  'TRP':'W',
  'ASP':'D',
  'GLU':'E',
  'ASN':'N', 
  'GLN':'Q',
  'LYS':'K',
  'ARG':'R',
  'HIS':'H',
  'ACE':'' ,
  'NH2':''
}


def filter_cycopep():


    # 遍历所有pdb文件
    sequence_lengths = []
    sequences = []
    for file in Path(work_path).glob('*.pdb'):
        if 'match' in file.name:

            # 打开pdb文件 
            with open(file) as f:
                lines = f.readlines()
            #print(file)
            # 找到SG原子所在行
            sg_atoms = []
            for line in lines:
                if line.startswith('HETATM') and line[12:16].strip() == 'SG':
                    sg_atoms.append(line)
            #print(sg_atoms)        
            if len(sg_atoms) == 2:
                # 提取SG坐标
                x1, y1, z1 = [float(sg_atoms[0][30:38]), float(sg_atoms[0][38:46]), float(sg_atoms[0][46:54])]
                x2, y2, z2 = [float(sg_atoms[1][30:38]), float(sg_atoms[1][38:46]), float(sg_atoms[1][46:54])]
                
                # 计算距离
                dist = ((x1-x2)**2 + (y1-y2)**2 + (z1-z2)**2)**0.5
                
                if dist < 2.05:
                    
                    # 提取序列                
                    sequence = ''
                    residue_ids = [] 
                    for line in lines:
                        if line.startswith('ATOM') or line.startswith('HETATM'):
                            res_id = line[22:27].strip()

                            res_key = res_id
                            if res_key not in residue_ids:
                                residue_ids.append(res_key)
                                res_name = line[17:20]
                                one_letter = residues_dict.get(res_name)
                                sequence += one_letter
                    sequences.append(sequence)
                    print(file, sequence,len(sequence))

                    # 写入文件
                    sequences.append(sequence)
                    sequence_lengths.append(len(sequence))

    plt.hist(sequence_lengths, bins=40,edgecolor='black')
    plt.xticks(range(0, 32, 2))
    plt.xlabel('Sequence length')  
    plt.ylabel('Count')
    plt.title('Sequence Length Distribution')
    plt.show()    
    return sequences

'''                
                with open('output.txt','a') as f:
                    f.write(file + '\n')
                    f.write(sequence + '\n')
'''
def colabfold(sequences):
    with open(protein_fasta) as f:
        protein_seq = ''
        for line in f:
            if not line.startswith('>'):
                protein_seq += line.strip()

    for peptide in sequences:
        if len(peptide) < 12:
            input_fasta_file = os.path.join(input_dir, peptide + '.fasta')
            if not os.path.exists(input_dir):
                os.makedirs(input_dir)
            with open(input_fasta_file, 'w') as f:
                f.write('>'+ peptide+'\n')
                f.write(protein_seq + ':\n')
                f.write( peptide)


    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
# python 解释器内核换成colabfold
    os.system('colabfold_batch  {}  {}  --amber --zip'.format(input_dir, output_dir))

if __name__ == '__main__':
    work_path  = '/mnt/sdc/lanwei/software/cPEPmatch/5w59/'
    protein_fasta = '/mnt/sdc/lanwei/software/cPEPmatch/5w59/fgf5.fasta'
    input_dir = '/mnt/sdc/lanwei/software/cPEPmatch/5w59/colab/input'
    output_dir = '/mnt/sdc/lanwei/software/cPEPmatch/5w59/colab/output'
    sequences = filter_cycopep()
    colabfold(sequences)