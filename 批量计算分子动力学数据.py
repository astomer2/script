import os
import csv
import pandas as pd
import matplotlib.pyplot as plt


def draw_RMSD_line():
    
    # 设置图形大小
    plt.figure(figsize=(10,8)) 

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
    plt.legend(loc = 'upper center',bbox_to_anchor=(0.5, 1.15), ncol = 5 )

    # 保存图像
    plt.savefig(RMSD_line,dpi = 600)

def plot_RMSD_box():
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
    plt.xticks(rotation=-45)
    plt.ylabel('RMSD (nm)')
    plt.grid(False)
    plt.savefig(RMSD_box,dpi = 600)

def take_energy():
    with open(Energy, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        
        writer.writerow(['Sequence', 'Energy'])
        
        for subdir in os.listdir(root_dir):
            filepath = os.path.join(root_dir, subdir, 'energy.dat')
            
            if os.path.exists(filepath):
                with open(filepath) as f:
                    for line in f:
                        if line.startswith('DELTA TOTAL'):
                            fields = line.split()
                            avg = fields[2]
                            std = fields[3]
                            
                            energy = avg + '±' + std
                            
                            writer.writerow([subdir, energy])
            
if __name__ == '__main__':
    root_dir = '/mnt/nas1/lanwei-125/IL8'
    RMSD_line = f'{root_dir}/RMSD_line.png'
    RMSD_box = f'{root_dir}/RMSD_box.png'
    Energy = f'{root_dir}/energy.csv'
    take_energy()
    draw_RMSD_line()
    plot_RMSD_box()
