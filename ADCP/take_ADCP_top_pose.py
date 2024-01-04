import os
import shutil


def confirm_best_pose_pdb():
    score_line = {}
    for root, dirs, files in os.walk(path):
        for dir in dirs:
            sub_path = os.path.join(root, dir)
            if 'log.txt' not in os.listdir(sub_path):
                continue
            log_file = open(sub_path + '/log.txt', 'r')
            if sum(1 for line in log_file) < 25:
                print(f'your docking job {dir} was failed,please cheak or rerun')
                continue
            log_file.seek(0)
            lines = log_file.readlines()
            for i in range(14, 25):
                if lines[i].split()[0] == '1':
                    score = lines[i].split()[1]
            
            score_line[dir] = score
    #pprint(score_line)

    min_scores  = {}
    for dirs, scores in score_line.items():
        sequence = (dirs.split("-")[0])
        # 如果sequence 不在min_score字典中，则将对应的分数写入字典
        # 如果sequence在字典中，则用新的分数和原来的分数做比较
        # 如果新的分数小于原来的分数，则为sequence对应的分数更新
        # 这样，可以保证sequence分数一直是最低的
        #在输出分数时，同时输出键的方向，则需要将输出的值转化为元组
        # 并且在比较分数时，从元组中提取数值，这要求在用float()时必须要指定值，而不是整个元组
        if sequence not in min_scores or float(scores) < float(min_scores[sequence][0]):
            min_scores[sequence] = (scores, dirs)
    #pprint(min_scores)
    return min_scores



def take_best_pose_pdb(min_scores):
    os.makedirs(out_dir, exist_ok=True)
    for peptide ,value in min_scores.items():
        score, folder = value
        sequence = peptide.upper()
        src_path = os.path.join(path +'/'+ folder , 'dock_ranked_1.pdb')

        if os.path.exists(src_path):
                dst_path = os.path.join(out_dir, sequence + '.pdb')
                shutil.copy(src_path, dst_path)
            #对PDB文件重排序
                with open(dst_path) as f:
                    lines = f.readlines()
                # 创建一个字典,键是残基编号,值是该残基的所有原子信息
                residues = {}
                for line in lines:
                    if line.startswith('ATOM'):
                        res_num = int(line[22:26])
                        if res_num not in residues:
                            residues[res_num] = []
                        residues[res_num].append(line)

                # 对每个残基进行排序
                for res_num, atom_lines in residues.items():
                    main_chain = []
                    side_chain = []
                    for line in atom_lines:
                        atom_name = line[12:16].strip()
                        if atom_name in ['N', 'CA', 'C', 'O']:
                            main_chain.append(line)
                        else:
                            side_chain.append(line)
                    
                    # 主链原子在前,侧链原子在后
                    sorted_lines = main_chain + side_chain
                    
                    # 重新编号
                all_lines = []
                    # 重新编号
                for res_num, lines in residues.items():
                    all_lines.extend(lines)
                    for i, line in enumerate(all_lines):
                        line = line[:6] + str(i+1).rjust(5) + line[11:]
                        all_lines[i] = line 

                # 把排序后的所有原子信息合并写入新PDB文件    
                with open(dst_path, 'w') as f:
                        f.writelines(all_lines)
                f.close()
                
            #对PDB文件进行修改
                with open(dst_path) as f:
                    content = f.readlines()
                    score_line = f'REMARK SCORE {score}\n'
                    content.insert(0, score_line)
                    #在修改行时，需要将新的行保存到变量中，这样在后续写入时才会保存修改
                    #用索引i来更新修改后的行到内容列表
                    for i, line in enumerate(content):
                        if line.startswith('ATOM'):

                            residue = line[17:20]
                            if residue == 'CYX':
                                line = line[:17] + 'CYS' + line[20:]
                            if residue == 'HIE':
                                line = line[:17] + 'HIS' + line[20:]

                            occupancy = line[54:60]
                            if occupancy != '1.00':
                                line = line[:56] + '1.00' + line[60:]

                            chain_id = line[21]  
                            if chain_id == 'A':
                                line = line[:21] + new_chain_id + line[22:]
                                content[i] = line
                                
                with open(dst_path, 'w') as f:
                    
                    f.writelines(content)

def merge_best_pdb():
    pdb_files = [f for f in os.listdir(out_dir) if f.endswith('.pdb')] 
    pdb_files.sort()

    with open(output_file, 'w') as f_out:
      for i, pdb_file in enumerate(pdb_files):
          with open(os.path.join(out_dir, pdb_file)) as f_in:
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
    path =  '/mnt/nas1/lanwei-125/FGF5/FGF5_disulfied_sequence_dock/'
    new_chain_id = 'E'
    out_dir = '/mnt/nas1/lanwei-125/FGF5/disf_ADCP_best_pose_pdb/'
    output_file = f'{out_dir}/MC1R-dock-pose.pdb'
    #md_file = '/mnt/nas1/lanwei-125/IL8/v2/'
    min_scores = confirm_best_pose_pdb()
    take_best_pose_pdb(min_scores)
    #merge_best_pdb()