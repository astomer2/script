import os
from re import T
import zipfile
from collections import defaultdict

def extract_zip_files(zip_folder, extracted_folder_base):
    """解压缩文件： extract_zip_files 函数遍历给定的文件夹 zip_folder 中的所有文件，
    如果文件的扩展名是 '.zip'，则将其解压到一个以文件名为基础的目标文件夹中。
    解压后的文件夹存储在 extracted_folder_base 文件夹中。"""
    for file_name in os.listdir(zip_folder):
        if file_name.endswith('.zip'):
            file_path = os.path.join(zip_folder, file_name)
            folder_name = os.path.splitext(file_name)[0]
            extracted_folder_path = os.path.join(extracted_folder_base, folder_name)
            os.makedirs(extracted_folder_path, exist_ok=True)
            
            with zipfile.ZipFile(file_path, 'r') as zip_ref:
                zip_ref.extractall(extracted_folder_path)

def find_sequences_and_rename(extracted_folder_base, cyclic_folder, pep_chain_id):
    """查找序列并重命名： find_sequences_and_rename 函数在 extracted_folder_base 文件夹及其子文件夹中查找具有特定命名模式的文件，
    提取相应的序列信息，根据一些条件（比如两个SG原子的距离），
    选择某个文件作为代表，并将其重命名为该序列的新名称，最后将这些文件移动到 cyclic_folder 文件夹中。"""
    sequences = defaultdict(lambda: (float('inf'), None))
    os.makedirs(cyclic_folder,exist_ok=True)
    for root, dirs, files in os.walk(extracted_folder_base):
        for file_name in files:
            if '_relaxed_rank_' in file_name and file_name.endswith('.pdb'):
                file_path = os.path.join(root, file_name)

                # Extract SG coordinates
                with open(file_path) as f:
                    lines = f.readlines()

                sg_atoms = [line for line in lines if line.startswith('ATOM') and line[21] == pep_chain_id and line[12:16].strip() == 'SG']

                if len(sg_atoms) == 2:
                    x1, y1, z1 = [float(sg_atoms[0][30:38]), float(sg_atoms[0][38:46]), float(sg_atoms[0][46:54])]
                    x2, y2, z2 = [float(sg_atoms[1][30:38]), float(sg_atoms[1][38:46]), float(sg_atoms[1][46:54])]
                    dist = ((x1 - x2) ** 2 + (y1 - y2) ** 2 + (z1 - z2) ** 2) ** 0.5

                    if dist < 2.05:
                        parts = file_name.split('_')
                        sequence = parts[0]
                        rank = int(parts[3])

                        if rank < sequences[sequence][0]:
                            sequences[sequence] = (rank, file_path)

    for seq, (rank, file_path) in sequences.items():
        outfile = os.path.join(cyclic_folder, f'{seq}.pdb')
        os.rename(file_path, outfile)

if __name__ == "__main__":
    sturcture_zip_path = '/mnt/nas1/lanwei-125/FGF5/FGF5-pos/new_pos/output/'
    extracted_folder_base = '/mnt/nas1/lanwei-125/FGF5/FGF5-pos/new_pos/MD/pos/'
    cyclic_folder = '/mnt/nas1/lanwei-125/FGF5/FGF5-pos/new_pos/MD/pos/'
    pep_chain_id = 'B'
    extract_zip_files(sturcture_zip_path, extracted_folder_base)
    find_sequences_and_rename(extracted_folder_base, cyclic_folder, pep_chain_id)
