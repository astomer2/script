import os
import csv
from tkinter import NO
import pandas as pd
import matplotlib.pyplot as plt 
from pathlib import Path
from numpy import around
import sys 
import argparse

def get_residue_num(pdb_path):
    """
    Counts the number of unique residues in a PDB file.

    Args:
        pdb_path (str): The path to the PDB file.

    Returns:
        int: The number of unique residues in the PDB file.
    """
    with open(pdb_path) as f:
        lines = f.readlines()

    residue_count = 0
    prev_res_id = 0
    peptide_id = 0 
    for line in lines:
        if line.startswith('ATOM'):
            peptide_id = line[22:26].strip()
            if peptide_id != prev_res_id:
                residue_count += 1
            prev_res_id = peptide_id
    return residue_count

def draw_RMSD_line(root_dir, RMSD_line):
    
    # 设置图形大小
    plt.figure(figsize=(15,8)) 

    for subdir in os.listdir(root_dir):
        if os.path.isdir(os.path.join(root_dir, subdir)):
            # 读取xvg文件
            filename = os.path.join(root_dir, subdir, 'complex_rmsd.xvg')
            
            # 判断xvg文件是否存在
            if not os.path.exists(filename):
                continue
            x = []
            y = []
            with open(filename,'r') as f:
                next(f) # 跳过第一行
                for line in f:
                    cols = line.split()
                    x.append(float(cols[0])/500) 
                    y.append(float(cols[1])/10)
            
            # 绘制折线图
            #plt.plot(x, y, '-', label=subdir) 
            #plt.plot(x, y, 'r-', label=subdir)       
            plt.plot(x, y, label=subdir, linewidth=0.5)

    # 设置轴标签        
    plt.xlabel('Time (ns)')
    plt.ylabel('RMSD (nm)') 

    # 设置图例
    plt.legend(loc = 'upper center',bbox_to_anchor=(0.5, 1.25), ncol = 5 )
    plt.tight_layout() 
    # 保存图像
    plt.savefig(RMSD_line,dpi = 600)

def plot_RMSD_box(root_dir, RMSD_box):
    plt.figure(figsize=(15,12))
    # 存储所有RMSD值的列表
    all_rmsd = []
    valid_folders = []

    # 遍历父文件夹下的所有子文件夹
    for sub_folder in os.listdir(root_dir):
        # 拼接每个子文件夹的路径
        file_path = os.path.join(root_dir, sub_folder, 'complex_rmsd.xvg')
        if not os.path.exists(file_path):
            continue
        # 读取RMSD值到列表
        rmsd = []
        with open(file_path) as f:
            for line in f:
                if line.startswith('#'):
                    continue
                cols = line.split()
                rmsd.append(float(cols[1])/10)          
        # 添加到总的RMSD列表中  
        all_rmsd.append(rmsd)
        valid_folders.append(sub_folder)
        
    # 将RMSD值转换为DataFrame
    df = pd.DataFrame(all_rmsd).T
    df.columns = valid_folders
    # 绘制箱线图
    plt.figure()
    df.boxplot(flierprops = {'marker':'o', 'markersize':1})
    plt.xticks(rotation=-90,size = 7)
    plt.ylabel('RMSD (nm)')
    plt.grid(False)
    plt.tight_layout() 
    plt.savefig(RMSD_box,dpi = 600)

def take_energy(Energy):

    with open(Energy, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        
        writer.writerow(['Sequence', 'avg', 'std'])
        
        for subdir in os.listdir(root_dir):
            filepath = os.path.join(root_dir, subdir, 'energy.dat')
            if os.path.exists(filepath):
                with open(filepath) as f:
                    for line in f:
                        if line.startswith('DELTA TOTAL'):
                            fields = line.split()
                            avg = fields[2]
                            std = fields[3]  
                            writer.writerow([subdir, avg, std])

def draw_energy_error_plot(Energy, energy_error_plot):
    energies = []  
    sequences = []
    with open(Energy) as f:
        reader = csv.reader(f)
        next(reader)
        for line in reader:
            sequence = line[0]
            avg = line[1]
            std = line[2]
            sequences.append(sequence)            
            energies.append((float(avg), float(std)))
    # 绘制误差线图 

    fig, ax = plt.subplots()

    ax.errorbar(sequences, [avg for avg, std in energies], 
                yerr=[std for avg, std in energies], fmt='o',  capsize= 5 )
    ax.set_xlabel('Sequence')  
    ax.set_ylabel('Energy (Kcal/mol)')
    plt.xticks(rotation=-90,size = 7)
    plt.tight_layout() 
    plt.savefig(energy_error_plot, dpi=600)

def residue_contribution(root_dir):
    root_dir = Path(root_dir)
    for subdir in root_dir.iterdir():
        if subdir.is_dir():
            if not os.path.exists(subdir / 'decomp-energy.dat'):
                continue
            else:
                file_path = subdir / 'decomp-energy.dat'
                data = pd.read_csv(file_path, skiprows=8, delimiter=',',header=None, engine='python')

                # 将第 0 列拆分为两个列
                a = data[0].str.split(expand=True, n=1)
                data['Residue ID'] = a[0]
                data['ids'] = a[1]

                # 将第 1 列拆分为三个列
                residue_info = data[1].str.split(expand=True, n=2)
                data['chain ID'] = residue_info[0]
                data['Amino Acid'] = residue_info[1]
                data['Residue_ID'] = residue_info[2]

                # 删除第 0 列和第 1 列
                data = data.drop([0, 1], axis=1)
                data["energy"] = data[5]

                # 提取所需列（以 Residue 列为 x 轴，Avg. 列为 y 轴）
                data_subset = data.iloc[get_residue_num(subdir/'protein.pdb'):]
                residue_ids = data_subset['Residue_ID']
                energy_contributions = data_subset['energy']

                # 绘制图表
                plt.plot(residue_ids, energy_contributions, label=subdir)
                plt.title('Energy Contributions of Amino Acid Residues')
                plt.xlabel('Peptide Residue ID')
                plt.ylabel('Energy Contribution (Kcal/mol)')
                plt.legend()
                plt.tight_layout() 
                plt.savefig(subdir / 'residue_contribution.png')
                plt.close()

def read_contacts(path,columns):
    nativate_contact = pd.read_csv(path, header=None)
    residue_info = nativate_contact[0].str.split(expand=True, n=4)

    residue_info.columns = residue_info.iloc[0]
    residue_info = residue_info[1:]
    filtered_df = residue_info[residue_info['TotalFrac'].astype(float) > 10]
    contact_res= (filtered_df.iloc[:,columns].astype(int)).tolist()
    
    return set(contact_res)

def contact(root_dir,ref_contact_data,contact_rate):
    root_dir = Path(root_dir)

    ref = Path(ref_contact_data)
    ref = set(x-225 for x in read_contacts(ref,1))

    rates = []
    for path in root_dir.iterdir():
        if path.is_dir():
            contact_dat = path/'contact_frac_byres.dat'
            if contact_dat.exists():
                contact = read_contacts(contact_dat,0)
                rate = len(contact & ref)/len(ref)
                a= path.stem , around(rate,4)
                rates.append(a)
    with open (contact_rate, 'w') as f:
        csv_writer = csv.writer(f)
        
        # 写入表头
        csv_writer.writerow(['Folder', 'Rate'])

        # 写入数据
        csv_writer.writerows(rates)

def analysis_MD_result(root_dir,ref_contact_data) -> None:

    RMSD_line = f'{root_dir}/RMSD_line.png'
    RMSD_box = f'{root_dir}/RMSD_box.png'
    Energy = f'{root_dir}/energy.csv'
    energy_error_plot = f'{root_dir}/energy_error_plot.png'
    contact_rate = f'{root_dir}/contact_rate.csv'


    take_energy(Energy)
    draw_RMSD_line(root_dir, RMSD_line)
    plot_RMSD_box(root_dir, RMSD_box)
    draw_energy_error_plot(Energy,energy_error_plot)
    residue_contribution(root_dir)
    if ref_contact_data != '':
        contact(root_dir, ref_contact_data, contact_rate) #this maybe need ad-hoc modification

if __name__ == '__main__':
    root_dir ='/mnt/nas1/lanwei-125/MC5R/MD/AFD/'
    ref_contact_data =''
    analysis_MD_result(root_dir,ref_contact_data)
    


# if __name__ == '__main__':
#     parse=argparse.ArgumentParser(description='MD_analysis')

#     parse.add_argument('-i','--root_dir',type=str,help='root_dir')
#     parse.add_argument('-r','--ref_contact_data',default='',type=str,help='ref_contact_data,this part may need ad-hoc modification in code')

#     args=parse.parse_args()

#     root_dir = args.root_dir
#     ref_contact_data = args.ref_contact_data

#     analysis_MD_result(root_dir,ref_contact_data)
    