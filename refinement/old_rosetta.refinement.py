import os
import shutil
import re
import pandas as pd
from pandas import DataFrame
#按照pdb文件名称新建文件夹，注意需要是蛋白-肽复合体
def read_pdb():
    pdb_name = []
    for file_name in os.listdir(dir_path):
        if file_name.endswith('.pdb'):
            pdb_name.append(file_name.split('.')[0])
            new_dir = os.path.join(dir_path, file_name.split('.')[0])
            if not os.path.exists(new_dir):
                os.makedirs(new_dir)
            pdb_file = os.path.join(dir_path, file_name)
            shutil.copy(pdb_file, new_dir)
    return pdb_name
    
#结构预处理，将残基命名，原子参数等设置正确，rosetta只能读取标注氨基酸命名，
#如果必须包含非标准氨基酸，则需要在flags文件中写明
def preprocess_pdb(pdb_name):
    for peptide in pdb_name:
        work_dir = os.path.join(dir_path, peptide)
        os.chdir(work_dir)
        #os.system("pdb4amber -i "+peptide+".pdb -o "+peptide+"-h.pdb --add-missing-atoms --reduce")
        pdb_name = os.path.join(dir_path, peptide, peptide+'.pdb')
        new_pdb = []
        with open(pdb_name) as f:
            for line in f:
                if line.startswith('ATOM'):
                    residue = line[17:20]
                    if residue == 'CYX':
                        line = line[:17] + 'CYS' + line[20:]
                    if residue == 'HIE':
                        line = line[:17] + 'HIS' + line[20:]
                    occupancy = line[54:60]
                    if occupancy != '1.00':
                        line = line[:56] + '1.00' + line[60:]
                new_pdb.append(line)

        with open(pdb_name, 'w') as f:
            f.write(''.join(new_pdb))


def make_flag(pdb_name):
    for peptide in pdb_name:
        work_dir = os.path.join(dir_path, peptide)
        os.chdir(work_dir)
        f=open('prepack.flags', 'w')
        f.write("-s "+ peptide + ".pdb \n")
        f.write("-flexpep_prepack \n")
        f.write("-nstruct 1 \n")
        f.write("-ex1 \n")
        f.write("-ex2aro \n")
        f.write("-mute core.io.database \n")
        f.write("-mute core.util.prof \n")
        f.write("-mute protocols.jd2.JobDistributor \n")
        f.write("-scorefile pack.score.sc \n")
        f.close()
        os.system("mpirun FlexPepDocking.mpi.linuxgccrelease @prepack.flags")

        f=open('refinement.flags', 'w')
        f.write("-s "+ peptide + "_0001.pdb \n")
        f.write("-native "+ peptide + ".pdb \n")
        f.write("-flexPepDocking:receptor_chain "+protein_chain+" \n")
        f.write("-flexPepDocking:peptide_chain "+peptide_chain+" \n")
        f.write("-pep_refine \n")
        f.write("-use_input_sc \n")
        f.write("-nstruct 100 \n")
        f.write("-scorefile refinement.score.sc \n")
        f.write("-ex1 \n")
        f.write("-ex2aro \n")
        f.write("-mute core.io.database \n")
        f.write("-mute core.util.prof \n")
        f.write("-mute protocols.jd2.JobDistributor \n")
        f.close()
        os.system("mpirun FlexPepDocking.mpi.linuxgccrelease @refinement.flags")

# 结构预处理
def take_candidate(pdb_name):
    if not os.path.exists(result_path): 
        os.makedirs(result_path)
    for peptide in pdb_name:
        score_file = os.path.join(dir_path, peptide, 'refinement.score.sc')       

        with open(score_file, 'r') as f:
            headers = None
            rows = []
            lines = f.readlines()
            for line in lines:  
                if line.startswith('SEQUENCE:'):
                    continue 
                elif line.startswith('SCORE:'):
                    if headers is None:
                        headers = line.split()[1:]
                    else:
                        rows.append(line.split()[1:])
                else:
                    rows.append(line.split())                            
            df = pd.DataFrame(rows, columns=headers)
            #将df中的reweighted_sc数据类型转为浮点
            df['reweighted_sc'] = df['reweighted_sc'].astype(float)
            #df的排序
            df = df.sort_values(by='reweighted_sc', ignore_index=True) 
            # 取出reweighted_sc最小值对应的行
            min_row = df.loc[df['reweighted_sc'].argmin()]  
        candidate_path = os.path.join(dir_path, peptide,min_row['description']+'.pdb')
        shutil.copy(candidate_path, result_path)
        result_score_file = os.path.join(result_path, 'candidate.txt')
        
        if not os.path.exists(result_score_file):
            with open(result_score_file, 'w') as f:
                f.write("name\tscore\n")
                f.write(f"{min_row['description']}\t{min_row['reweighted_sc']}\n")
        else:
            with open(result_score_file, 'a') as f:
                f.write(f"{min_row['description']}\t{min_row['reweighted_sc']}\n")

# 结构细化
def process_refine_pdb():
    for file_name in os.listdir(result_path):
        if file_name.endswith('.pdb'):
            process_pdb_path = os.path.join(result_path, file_name)
            candidate_pdb_path = os.path.join(result_path, file_name.split('_')[0]+'.pdb')
            candidate_pdb = []
            with open(process_pdb_path) as f:
                lines = f.readlines()
                for i, line in enumerate(lines):
                    if line.startswith('ATOM'):
                        candidate_pdb.append(line)
                    if line.startswith('TER'):
                        candidate_pdb.append(line)
            with open(candidate_pdb_path, 'w') as f:
                f.write(''.join(candidate_pdb))


if __name__ == '__main__':
    dir_path= '/mnt/sdc/lanwei/IL8/test/'
    result_path = f'{dir_path}/refinement/'
    protein_chain = "A"
    peptide_chain = "C"
    pdb_name = read_pdb()
    preprocess_pdb(pdb_name)
    make_flag(pdb_name)
    take_candidate(pdb_name)
    process_refine_pdb()