import os
import shutil

input_file = r"C:\Users\123\Desktop\jobwork\subject\ADIPOR1\对接结果验证\hpep"
result_dir = r'C:\Users\123\Desktop\jobwork\subject\ADIPOR1\对接结果验证\hpep\hpep-result'

os.makedirs(result_dir, exist_ok=True)

for root, _, files in os.walk(input_file):

    # 解析每个文件夹的hpepdock_all.out
    peptide_info = {}
    if 'hpepdock_all.out' in files:
        with open(os.path.join(root, 'hpepdock_all.out')) as f:
            for line in f:
                data = line.split()
                peptide = data[4]
                cluster = data[0]
                number = data[1]
                score = float(data[3])
                if peptide not in peptide_info:
                    peptide_info[peptide] = []
                peptide_info[peptide].append((cluster, number, score))

        # 在该文件夹下提取最佳模型            
    for peptide, info in peptide_info.items():
        cluster, number, _ = info[0]  # 取第一个cluster和number
        pdb_file_path = os.path.join(result_dir, f'{peptide}.pdb')

        # 读取 hpepdock_all.pdb 文件
        pdb_data = []
        pdb_section = False
        if 'hpepdock_all.pdb' in files:
            with open(os.path.join(root, 'hpepdock_all.pdb')) as pdb_file:
                for line in pdb_file:
                    if line.startswith('REMARK Cluster:'):
                        # 检查当前结构的cluster和number是否匹配
                        _, file_cluster, file_number = line.split()[-3:]
                        if cluster == file_cluster and number == file_number:
                            pdb_section = True
                        else:
                            pdb_section = False
                    elif pdb_section:
                        pdb_data.append(line)

            # 保存提取的pdb结构
            with open(pdb_file_path, 'w') as output_pdb:
                output_pdb.writelines(pdb_data)

print("PDB文件提取完成")