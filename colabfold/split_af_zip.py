import os
import zipfile
import shutil
import sys
sys.path.append(os.path.abspath('.'))
from MD.split_pdb_file import split_pdb_run

def extract_and_move_pdb_files(zip_folder, output_folder):
    os.makedirs(output_folder, exist_ok=True)

    for file_name in os.listdir(zip_folder):
        if file_name.endswith('.zip'):
            file_path = os.path.join(zip_folder, file_name)

            with zipfile.ZipFile(file_path, 'r') as zip_ref:
                # 创建一个目录用于解压缩文件
                extract_folder = os.path.join(zip_folder, 'extracted_files')
                os.makedirs(extract_folder, exist_ok=True)
                zip_ref.extractall(extract_folder)

                for root, dirs, files in os.walk(extract_folder):
                    for file in files:
                        if '_relaxed_rank_001' in file and file.endswith('.pdb'):
                            print(file)
                            name_parts = file.split('_')
                            pdb_name = name_parts[0] + '.pdb'
                            shutil.move(os.path.join(root, file), os.path.join(output_folder, pdb_name))

                # 删除临时的解压缩目录
                shutil.rmtree(extract_folder)

def organize_pdb_files_by_folder(input_folder):
    for pdb_name in os.listdir(input_folder):
        if pdb_name.endswith('.pdb'):
            print(pdb_name)
            pdb_name_list = os.path.splitext(pdb_name)[0]
            folder_path = os.path.join(input_folder, pdb_name_list)

            os.umask(0)
            # Create the folder if it doesn't exist
            os.makedirs(folder_path, exist_ok=True)

            # Move the .pdb file into its corresponding folder
            shutil.copy(os.path.join(input_folder, pdb_name), os.path.join(folder_path, pdb_name))

if __name__ == "__main__":
    sturcture_zip_path = '/mnt/sdc/lanwei/TGF/CF_relax/output/'
    top_pose_path = '/mnt/sdc/lanwei/TGF/CF_relax/extracted_files/'

    extract_and_move_pdb_files(sturcture_zip_path, top_pose_path)
    organize_pdb_files_by_folder(top_pose_path)
    split_pdb_run(top_pose_path)