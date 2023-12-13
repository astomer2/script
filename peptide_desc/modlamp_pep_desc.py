import os
import csv
import numpy as np
import pandas as pd
from modlamp.descriptors import GlobalDescriptor,PeptideDescriptor
from collections import OrderedDict

'''def calculate_descriptors(input_file, output_file):
     从输入文件中读取多肽序列
    with open(input_file, 'r') as file:
        sequences = file.readlines()
    unique_sequences = list(OrderedDict.fromkeys(sequences))
'''


#计算氨基酸二级结构倾向 https://link.springer.com/article/10.1007/s008940050118
#20个氨基酸二级结构倾向权重
aa_propensity = {
    'A': [1.41, 0.82, 0.72],
    'R': [1.21, 0.90, 0.84],
    'N': [0.76, 1.34, 0.48],
    'D': [0.99, 1.24, 0.39],
    'C': [0.66, 0.54, 1.40],
    'Q': [1.27, 0.84, 0.98],
    'E': [1.59, 1.01, 0.52],
    'G': [0.43, 1.77, 0.58],
    'H': [1.05, 0.81, 0.80],
    'I': [1.09, 0.47, 1.67],
    'L': [1.34, 0.57, 1.22],
    'K': [1.23, 1.07, 0.69],
    'M': [1.30, 0.52, 1.14],
    'F': [1.16, 0.59, 1.33],
    'P': [0.34, 1.32, 0.31],
    'S': [0.57, 1.22, 0.96],
    'T': [0.76, 0.96, 1.17],
    'W': [1.02, 0.65, 1.35],
    'Y': [0.74, 0.76, 1.45],
    'V': [0.90, 0.41, 1.87]
}

#计算多肽二级结构倾向
def get_sec_propensity(seq):
    alpha = 0
    turn = 0
    beta = 0
    for aa in seq:
        alpha += aa_propensity[aa][0]
        turn += aa_propensity[aa][1] 
        beta += aa_propensity[aa][2]
    
    length = len(seq)
    return [round(alpha/length, 3), round(turn/length, 3), round(beta/length, 3)]


#20个氨基酸logP(正辛醇/水分配系数)和logD（正辛醇/缓冲溶液分配系数）倾向权重

aa_dict = {'A':(-0.27,-0.27), 
           'R':(-0.79,-1.65), 
           'N':(-0.98,-0.98), 
           'D':(-0.28,-2.06), 
           'C':(0.83,0.82), 
           'Q':(-1.00,-1.00),
           'E':(-0.34,-2.19), 
           'G':(-0.22,-0.22), 
           'H':(-0.31,-0.44),
           'I':(0.70,0.69), 
           'L':(0.80,0.80), 
           'K':(0.17,-2.27), 
           'M':(0.51,0.51), 
           'F':(1.16,1.16), 
           'P':(0.15,0.15),
           'S':(-0.45,-0.45), 
           'T':(-0.26,-0.26), 
           'W':(1.46,1.46),
           'Y':(0.55,0.55), 
           'V':(0.32,0.32)}

Bp = -1.19 
Up = -3.25
BD = -1.18
UD = -3.25

#计算多肽logP和logD值
def calc_logP_logD(sequence):
    logP = 0
    logD = 0 
    if sequence[:3] == 'AC-' and sequence[-3:] == '-NH':
        is_protected = True
    else:
        is_protected = False

    for aa in sequence:
        logP += aa_dict[aa][0]
        logD += aa_dict[aa][1]
    
    if is_protected:
        logP += Bp
        logD += BD
    else:
        logP += Up
        logD += UD
    
    return logP, logD

#计算多肽的渗皮性能 https://analyticalsciencejournals.onlinelibrary.wiley.com/doi/full/10.1002/jat.4435

def calc_logKp_Kp(logD, MW):
    logKp = 0.65 * logD - 0.0085 * MW - 2.55
    Kp = 10**logKp

    return logKp, Kp

# 氨基酸交互矩阵 https://www.sciencedirect.com/science/article/pii/S1476927109001091
aa_matrix = np.array([[100, 100, 48.6, 41.8, 100, 100, 100, 100, 126, 130, 150, 58.8, 100, 36, 54.4, 100, 141, 168, 100, 87],
                      [28.2, 100, -78, 29.3, 100, 161, 100, 32.1, 100, 100, 100, 100, 255, 100, 100, 100, -62.1, 100, 100, 100],
                      [100, 100, 100, 100, 142, -12, 24.7, 55.8, 100, 134, 25.1, 52.4, 100, 165, 100, 143, 48, 140, -99, 100],
                      [100, 26.7, 100, 110, 100, 178, 100, 100, 100, 64.4, 102, 143, 62.6, 54.5, 142, 83.1, 115, 100, 100, 63.5],
                      [100, 168, 100, 100, 3.45, 100, 248, 87.2, 153, 16.7, 100, 100, 3.45, 168, 100, 144, 100, 26.3, 197, 100],
                      [100, 100, 100, 174, 107, 124, 100, 150, 66.5, 105, 51.9, 138, 100, -15, 61.3, 138, 100, 2.2, 151, 67.3],
                      [100, 402, 93, 100, 177, 100, 100, -59, 24.5, 33.4, 100, 36, 189, 100, 28.6, 184, 100, 100, 100, 100],
                      [100, -43, 35.3, 100, -37, 136, 244, 100, 100, 100, 100, 100, 100, 100, 100, 100, 158, 100, 100, 153],
                      [100, 29, 100, 100, 105, 100, 100, 100, 141, 47.5, 100, 165, 58.3, 100, 137, 57.5, 100, 132, 100, 121],
                      [100, 56.3, 100, 132, 100, 63.7, 100, 100, 100, 133, 32.2, 100, 100, 54.2, 100, 100, 100, 100, 100, 173],
                      [164, 100, 226, 23.4, 100, 153, 100, 100, 100, 100, -209, 219, 100, 100, 100, 46.6, -37.2, 37.1, 100, 100],
                      [100, -76, 100, 102, 100, 100, 100, 100, 3.49, 100, 164, -2.25, 100, 116, 166, 100, 100, 161, 100, 88.1],
                      [154, 159, 100, 100, 100, 61, 100, 100, 100, 116, 100, 100, 139, 168, 34, 100, 52.3, 88.6, 100, 100],
                      [100, -8.5, 147, 156, 42.1, 132, 188, 60.1, 100, 95.1, 100, 100, 100, 253, 15.8, 63.7, 100, 43.6, -40.6, 149],
                      [62.1, 100, 140, 66.3, 100, 100, 25.9, 16.3, 137, 100, 100, 100, 159, 47.5, 144, 67.6, 100, 93.6, 100, 100],
                      [68.8, 100, 159, 68.6, 149, 100, 88.5, 100, 100, 136, 100, 100, 100, 100, 100, 44.6, 100, 55.6, 100, 100],
                      [129, 100, 100, 106, 100, 100, 100, 100, 34.1, 100, 100, 100, 62.8, 100, 100, 100, 100, 100, -26.8, 100],
                      [100, 173, 100, 100, 100, 73.4, 24.5, 100, 174, 158, 100, 100, 43.6, 100, 53.3, 100, 87.3, 146, 100, 41.2],
                      [100, -153, 100, 25, 100, 151, 100, 211, 100, -69, 306, 100, -84, 100, 158, 25.5, 100, 100, 100, -23],
                      [100, 100, 151, 140, 100, 11.7, 191, 100, 100, 100, 100, 88.1, 197, 212, 100, 40.1, 100, 30.9, 100, 21.6]])

# 计算蛋白质熔解温度因子
def calculate_TI(sequences):
    amino_acids = 'ACDEFGHIKLMNPQRSTVWY'
    L = len(sequences)
    TI = 0
    for i in range(L-1):
        aa1 = sequences[i]
        aa2 = sequences[i+1]
        idx1 = amino_acids.index(aa1)
        idx2 = amino_acids.index(aa2)
        TI += aa_matrix[idx1, idx2]
    TI = ((100/L) * TI - 9372)/398
    return TI





def calculate_descriptors(input_file, output_file, columns_index, sheet_name):

    sequences = []
    
    # 判断文件类型
    if input_file.endswith('.csv'):
        # 读取CSV
        with open(input_file, 'r') as file:
            reader = csv.reader(file)
            for row in reader:
                sequences.append(row[columns_index])
        
    elif input_file.endswith('.txt'):
        # 读取TXT
        with open(input_file, 'r') as file:
            for line in file:
                sequences.append(line.strip())
    elif input_file.endswith('.xlsx'):
        # 读取XLSX
        df = pd.read_excel(input_file, sheet_name=sheet_name, header=0)
        sequences = df.iloc[:, columns_index].tolist()
    else:
        print("Unsupported file format!")
        return
    
    unique_sequences = list(OrderedDict.fromkeys(sequences))

    # 去除每个序列末尾的换行符
    #unique_sequences = [seq.strip() for seq in unique_sequences[1:]]

    # 创建CSV文件并写入表头
    with open(output_file, 'w', newline='') as file:
        writer = csv.writer(file)
        header = ['Sequence', 
                  'Length', 
                  'MW', 
                  'Charge', 
                  #'ChargeDensity', 
                  'pI', 
                  #'InstabilityInd', 
                  #'Aromaticity', 
                  #'AliphaticInd', 
                  #'BomanInd', 
                  #'HydrophRatio', 
                  'gravy',
                  'alpha-helix propensity',
                  #'side chains into lipid',
                  'Alipha-helix propensity',
                  'Reverse turn',
                  'Beta-sheet propensity',
                  'LogP',
                  'LogD',
                  'logKp',
                  'Kp',
                  'TI'


                  ]
        writer.writerow(header)

        # 逐个计算描述符并将结果写入CSV文件
        for unique_sequence in unique_sequences:
            desc = GlobalDescriptor(unique_sequence)
            desc.calculate_all(amide=True)

            pepdesc = PeptideDescriptor(unique_sequence)
            pepdesc.load_scale('gravy')
            pepdesc.calculate_global(append=True)

            alpha = PeptideDescriptor(unique_sequence)
            alpha.load_scale('levitt_alpha')
            alpha.calculate_global(append=True)

            #side_chains = PeptideDescriptor(unique_sequence)
            #side_chains.load_scale('Ez')
            #side_chains.calculate_global(append=True)

            sec_propensity = get_sec_propensity(unique_sequence)

            log_values = calc_logP_logD(unique_sequence)

            Skin_penetration_rate = calc_logKp_Kp(log_values[0], desc.descriptor[0][1]) 

            TM_index = calculate_TI(unique_sequence)

            row = [
                unique_sequence,
                desc.descriptor[0][0],       # Length
                desc.descriptor[0][1],       # MW
                desc.descriptor[0][2],       # Charge
                #desc.descriptor[0][3],       # ChargeDensity
                desc.descriptor[0][4],       # pI
                #desc.descriptor[0][5],       # InstabilityInd
                #desc.descriptor[0][6],       # Aromaticity
                #desc.descriptor[0][7],       # AliphaticInd
                #desc.descriptor[0][8],       # BomanInd
                #desc.descriptor[0][9],       # HydrophRatio
                pepdesc.descriptor[0][0],   # gravy
                alpha.descriptor[0][0],       # levitt_alpha
                #side_chains.descriptor[0][0], # side chains into lipid
                sec_propensity[0],       # Alipha-helix propensity
                sec_propensity[1],       # Reverse turn
                sec_propensity[2],        # Beta-sheet propensity
                log_values[0],        # LogP
                log_values[1],        # LogD
                Skin_penetration_rate[0],        # logKp
                Skin_penetration_rate[1],        # Kp
                TM_index       # TI,不用索引，直接输出

                
            ]
            writer.writerow(row)

    print(f"描述符已计算并保存到文件：{output_file}")


# 指定输入和输出文件的路径
if __name__ == "__main__":
    input_file = '/mnt/nas1/lanwei-125/TGFbR2/TGFBR2-adcp-result-sort.xlsx'
    output_file ='/mnt/nas1/lanwei-125/TGFbR2/TGFBR2-ph.csv'
    columns_index = 0
    sheet_name = 'phsical_properties'
    # 计算描述符并保存结果到CSV文件
    calculate_descriptors(input_file, output_file, columns_index, sheet_name)
