import os
import shutil

def collect_pdbfile(target_path, output_folder):
    # 移动当前目录到目标目录
    target_directory = target_path
    os.chdir(target_directory)

    # 获取当前目录下所有文件夹
    current_directory = os.getcwd()
    subdirectories = [folder for folder in os.listdir(current_directory) if os.path.isdir(folder)]

    # 遍历每个文件夹
    for folder in subdirectories:
        folder_path = os.path.join(current_directory, folder)
        os.chdir(folder_path)  # 进入文件夹

        # 检查文件夹中是否包含 protein.pdb 和 lig.pdb
        if "protein.pdb" in os.listdir(folder_path) and "lig.pdb" in os.listdir(folder_path):
            
            # 复制 protein.pdb
            protein_src = os.path.join(folder_path, "protein.pdb")
            protein_dst = os.path.join(output_folder, folder + "_protein.pdb")
            shutil.copy(protein_src, protein_dst)

            # 复制 peptide.pdb
            peptide_src = os.path.join(folder_path, "lig.pdb")
            peptide_dst = os.path.join(output_folder, folder + "_peptide.pdb")
            shutil.copy(peptide_src, peptide_dst)
            
            print(f"{folder} copied")

        os.chdir(current_directory)  # 返回上级目录
    print("finished")
        
if __name__ == "__main__":
    output_folder = "/mnt/nas/yuliu/repos/Dock2MD/pdb_files"
    target_paths = []

    target_paths.append("/mnt/sdc/hongbw/9.27/peptide-later30/2")
    target_paths.append("/mnt/sdc/hongbw/9.27/peptide-later30/1")
    target_paths.append("/mnt/sdc/hongbw/9.27/peptide-later30/4")
    target_paths.append("/mnt/sdc/hongbw/9.27/peptide-later30/3")
    target_paths.append("/mnt/sdc/hongbw/9.20/peptide/2")
    target_paths.append("/mnt/sdc/hongbw/9.20/peptide/1")
    target_paths.append("/mnt/sdc/hongbw/9.20/peptide/3")
    target_paths.append("/mnt/sdc/hongbw/23.8.16/MD")
    target_paths.append("/mnt/sdc/hongbw/9.13/peptide")
    target_paths.append("/mnt/sdc/hongbw/23.8.3/dock2/MD")
    target_paths.append("/mnt/sdc/hongbw/8.29/MD-zhaunli")

    '''
    target_paths.append("/mnt/nas1/lanwei-125/FGF5")
    target_paths.append("/mnt/nas1/lanwei-125/FGF5/FGF5-pos/new_pos/MD/pos")
    target_paths.append("/mnt/nas1/lanwei-125/FGF5/FGF5-pos/cpepmathc")
    target_paths.append("/mnt/nas1/lanwei-125/FGF5/Dis_to_cov/CC_to_SS/MD")
    target_paths.append("/mnt/nas1/lanwei-125/FGF5/Dis_to_cov/CC_to_KD")
    target_paths.append("/mnt/nas1/lanwei-125/FGF5/Dis_to_cov/CC_ro_AM/MD")
    '''
    for target_path in target_paths:
        collect_pdbfile(target_path, output_folder)
    
        
