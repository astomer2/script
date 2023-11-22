import os
import shutil

#将最优结构提取出来
def copy_ranked_pdbs(pdb_dir, top_pdb_dir):
    if not os.path.exists(top_pdb_dir):
        os.makedirs(top_pdb_dir)
        
    for root, dirs, files in os.walk(pdb_dir):
        for fname in files:
            if 'rank_001' in fname and fname.endswith('.pdb'):
                pdb_path = os.path.join(root, fname)
                folder_name = (os.path.basename(root)).upper()
                new_name = f'{folder_name}_ranked_1.pdb'
                top_pdb_path = os.path.join(top_pdb_dir, new_name)
                shutil.copy(pdb_path, top_pdb_path)
                #print(f'Copied {pdb_path} to {top_pdb_path}')
                

#合并最优结构
def combine_ranked_pdbs(top_pdb_dir, output_file):

    pdb_files = [f for f in os.listdir(top_pdb_dir) if f.endswith('.pdb')] 
    pdb_files.sort()

    with open(output_file, 'w') as f_out:
      for i, pdb_file in enumerate(pdb_files):
          with open(os.path.join(top_pdb_dir, pdb_file)) as f_in:
            f_out.write(f'\nMODEL {i+1}\n')
            
            lines = []
            for line in f_in:
                if 'MODEL' and 'END' and 'ENDMDL' not in line:
                    lines.append(line)
            f_out.write(''.join(lines))

            f_out.write(f_in.read())
            f_out.write('END\n')
            f_out.write('ENDMDL\n\n')

if __name__ == '__main__':
    pdb_dir = r'C:\Users\123\Desktop\jobwork\subject\IL-8\结构\colab'  
    top_pdb_dir = r'C:\Users\123\Desktop\jobwork\subject\IL-8\结构\colab\ranked'
    output_file = f'{top_pdb_dir}\\all_ranked_pdbs.pdb' 
    copy_ranked_pdbs(pdb_dir, top_pdb_dir)
    combine_ranked_pdbs(top_pdb_dir, output_file)





