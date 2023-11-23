import os
import zipfile
import shutil

def extract_and_move_pdb_files(zip_folder, output_folder):
    os.makedirs(output_folder, exist_ok=True)

    for file_name in os.listdir(zip_folder):
        if file_name.endswith('.zip'):
            file_path = os.path.join(zip_folder, file_name)

            with zipfile.ZipFile(file_path, 'r') as zip_ref:
                for f in zip_ref.namelist():
                    if '_relaxed_rank_001' in f and f.endswith('.pdb'):
                        print(f)
                        name_parts = f.split('_')
                        pdb_name = '_'.join(name_parts[:2]) + '.pdb'
                        zip_ref.extract(f, zip_folder)  
                        shutil.move(os.path.join(zip_folder, f), os.path.join(output_folder, pdb_name))

def organize_pdb_files_by_folder(input_folder):
    for pdb_name in os.listdir(input_folder):
        if pdb_name.endswith('.pdb'):
            print(pdb_name)
            pdb_name_list = os.path.splitext(pdb_name)[0]
            folder_path = os.path.join(input_folder, pdb_name_list)

            # Create the folder if it doesn't exist
            os.makedirs(folder_path, exist_ok=True, mode=0o777)

            # Move the .pdb file into its corresponding folder
            shutil.move(os.path.join(input_folder, pdb_name), os.path.join(folder_path, pdb_name))

if __name__ == "__main__":
    sturcture_zip_path = '/mnt/nas1/lanwei-125/IL8/v4/AF/'
    top_pose_path = '/mnt/nas1/lanwei-125/IL8/v4/MD/'

    extract_and_move_pdb_files(sturcture_zip_path, top_pose_path)
    organize_pdb_files_by_folder(top_pose_path)
