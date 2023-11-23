import os
import shutil
import re
import pandas as pd
from pandas import DataFrame
from icecream import ic
'''
#按照pdb文件名称新建文件夹，注意需要是蛋白-肽复合体
def read_pdb():
    pdb_names = []
    for file_name in os.listdir(dir_path):
        if file_name.endswith('.pdb'):
            pdb_names.append(file_name.split('.')[0])
            new_dir = os.path.join(dir_path, file_name.split('.')[0])
            if not os.path.exists(new_dir):
                os.makedirs(new_dir)
            pdb_file = os.path.join(dir_path, file_name)
            shutil.copy(pdb_file, new_dir)
    return pdb_names
'''    
# 此脚本将与聚类结果偶联，如果没做或者没有聚类文件夹,则启用上面被注释的代码，记得修改下面的路径命名

#读取pdb文件所在的工作路径

def read_pdb():
    cluster_dir_names = []
    pdb_names = []
    for dir_name in os.listdir(dir_path):
        if os.path.isdir(os.path.join(dir_path, dir_name)) and dir_name.split("_")[0] == "cluster":
            cluster_dir_names.append(dir_name)

    for dir_name in cluster_dir_names:
        pdb_path = os.path.join(dir_path, dir_name)
        for file in os.listdir(pdb_path):
            if file.split(".")[1] == "pdb":
                pdb_names.append(file.split(".")[0])
    return cluster_dir_names, pdb_names

#结构预处理，将残基命名，原子参数等设置正确，为多肽添加蛋白质，蛋白质添加在多肽之前。rosetta只能读取标注氨基酸命名，
#如果必须包含非标准氨基酸，则需要在flags文件中写明
def preprocess_pdb(protein_pdb, cluster_dir_names):
    protein_lines = []
    with open(protein_pdb) as f:
        protein_liners = f.readlines()
        for protein_line in protein_liners:
            if protein_line.startswith("ATOM"):
                protein_lines.append(protein_line)
            if protein_line.startswith("TER"):
                protein_lines.append(protein_line)

    for cluster_dir_name in cluster_dir_names:
        pdb_path = os.path.join(dir_path, cluster_dir_name)
        
        for file in os.listdir(pdb_path):
            pdb_file = os.path.join(pdb_path, file)
            parts = file.split("-")
            if len(parts) >= 2 and parts[1] == "-ref.pdb":
                continue

            #在写if/elif时条件不要形成交集
            elif file.endswith('.pdb') and '-ref' not in file:
                # 读取对接结果的pdb文件
                peptide_lines = []
                with open(pdb_file) as f:
                    lines = f.readlines()
                    for line in lines:
                        if line.startswith("ATOM"):
                            peptide_lines.append(line)

                new_pdb_path = os.path.join(pdb_path, file.split(".")[0]+'-ref.pdb')
                with open(new_pdb_path, 'w') as f:
                    f.write(''.join(protein_lines))                  
                    if not protein_lines[-1].startswith('TER'):
                        f.write("\nTER\n")
                    f.write(''.join(peptide_lines))
                    f.write("\nTER\n")
                    f.write("END\n")

# 使用rosetta_flex_pepdock进行结构细化处理
# 首先预处理pdb文件，然后进行局部细化，输出打分文件refinement.score.sc
def make_flag(cluster_dir_names):
    for cluster_dir_name in cluster_dir_names:
        work_dir = os.path.join(dir_path, cluster_dir_name)
        os.chdir(work_dir)
        ic(work_dir)

        for file in os.listdir(work_dir):
            ic(file)
            parts = file.split("-")
            if len(parts) >= 2 and parts[1] == "ref.pdb":
                ref_pdb = file
                peptide = ref_pdb.split(".")[0]

                f=open('prepack.flags', 'w')
                f.write("-s "+ peptide + ".pdb \n")
                f.write("-flexpep_prepack \n")
                f.write("-nstruct 1 \n")
                f.write("-ex1 \n")
                f.write("-ex2aro \n")
                f.write("-mute core.io.database \n")
                f.write("-mute core.util.prof \n")
                f.write("-mute protocols.jd2.JobDistributor \n")
                f.write("-scorefile pack.score.sc ")
                f.close()
                os.system("mpirun -np 30 FlexPepDocking.mpi.linuxgccrelease @prepack.flags")

                f=open('refinement.flags', 'w')
                f.write("-s "+ peptide + "_0001.pdb \n")
                #f.write("-native "+ peptide + ".pdb \n")
                f.write("-flexPepDocking:receptor_chain "+protein_chain+" \n")
                f.write("-flexPepDocking:peptide_chain "+peptide_chain+" \n")
                f.write("-pep_refine \n")
                f.write("-use_input_sc \n")
                f.write("-nstruct 250 \n")
                f.write("-scorefile refinement.score.sc \n")
                f.write("-ex1 \n")
                f.write("-ex2aro \n")
                f.write("-mute core.io.database \n")
                f.write("-mute core.util.prof \n")
                f.write("-mute protocols.jd2.JobDistributor ")
                f.close()
                os.system("mpirun -np 30 FlexPepDocking.mpi.linuxgccrelease @refinement.flags")

def take_candidate(cluster_dir_names, pdb_names):
    if not os.path.exists(result_path): 
        os.makedirs(result_path)
    for cluster_dir_name in cluster_dir_names:
        for peptide in pdb_names:
            score_file = os.path.join(dir_path, cluster_dir_name, 'refinement.score.sc')       

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

            candidate_path = os.path.join(dir_path, cluster_dir_name, min_row['description']+'.pdb')
            shutil.copy(candidate_path, result_path)

        result_score_file = os.path.join(result_path, 'candidate.txt')
        if not os.path.exists(result_score_file):
            with open(result_score_file, 'w') as f:
                f.write("name\tscore\n")
                f.write(f"{min_row['description']}\t{min_row['reweighted_sc']}\n")
        else:
            with open(result_score_file, 'a') as f:
                f.write(f"{min_row['description']}\t{min_row['reweighted_sc']}\n")

# 局部细化后处理
def process_refine_pdb():
    for file_name in os.listdir(result_path):
        if file_name.endswith('.pdb'):
            processed_pdb_path = os.path.join(result_path, file_name)
            candidate_pdb_path = os.path.join(result_path, file_name.split('-')[0]+'.pdb')
            candidate_pdb = []
            with open(processed_pdb_path) as f:
                lines = f.readlines()
                for i, line in enumerate(lines):
                    if line.startswith('ATOM'):
                        candidate_pdb.append(line)
                    if line.startswith('TER'):
                        candidate_pdb.append(line)
                    if line.startswith('END'):
                        candidate_pdb.append(line)
            with open(candidate_pdb_path, 'w') as f:
                f.write(''.join(candidate_pdb))


if __name__ == '__main__':
    dir_path= '/mnt/sdc/lanwei/MC1R/final_best_pose_pdb/'
    result_path = f'{dir_path}/refinement/'
    protein_pdb = '/mnt/sdc/lanwei/MC1R/MC1R.pdb'
    protein_chain = "A"
    peptide_chain = "E"
    cluster_dir_names, pdb_names = read_pdb()
    preprocess_pdb(protein_pdb, cluster_dir_names)
    make_flag(cluster_dir_names)
    take_candidate(cluster_dir_names, pdb_names)
    process_refine_pdb()