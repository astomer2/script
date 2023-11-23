import os
import shutil

def copy_ranked_pdbs(src_dir, dest_dir):
  if not os.path.exists(dest_dir):
    os.makedirs(dest_dir)
  
  for root, dirs, files in os.walk(src_dir):
    for fname in files:
      if 'ranked_1' in fname and fname.endswith('.pdb'):
        src_path = os.path.join(root, fname)
        folder_name = (os.path.basename(root)).upper()
        new_name = f'{folder_name}_ranked_1.pdb'
        dest_path = os.path.join(dest_dir, new_name)
        shutil.copy(src_path, dest_path)
        print(f'Copied {src_path} to {dest_path}')

if __name__ == '__main__':
  src_dir = '/mnt/sdc/lanwei/IL8/colabfold/output/'
  dest_dir = '/mnt/sdc/lanwei/IL8/colabfold/output/ranked'
  copy_ranked_pdbs(src_dir, dest_dir)