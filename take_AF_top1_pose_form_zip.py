import os
import zipfile
import shutil

original_folder = '/mnt/sdc/lanwei/ADIPOR/colabfold/output/'
new_folder = '/mnt/sdc/lanwei/ADIPOR/colabfold/top-pdb'

if not os.path.exists(new_folder):
    os.makedirs(new_folder)

for file_name in os.listdir(original_folder):
    if file_name.endswith('.zip'):
        file_path = os.path.join(original_folder, file_name)
        
        with zipfile.ZipFile(file_path, 'r') as zip_ref:
            for f in zip_ref.namelist():
                if 'rank_001' and 'pdb' in f:
                    name_parts = f.split('_')
                    pdb_name = name_parts[0] + '.pdb'
                    zip_ref.extract(f, original_folder)  
                    shutil.move(os.path.join(original_folder, f), os.path.join(new_folder, pdb_name))

os.chdir(original_folder)
os.system('zip -r top-pdb.zip .')