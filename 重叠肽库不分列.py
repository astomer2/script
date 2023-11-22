#重叠肽库，不分列版本

import csv
from collections import OrderedDict
from pprint import pprint

# 肽链名称和序列
peptide_sequences = [
    {"name": "HBD3", 
     "sequence": "GIINTLQKYYCRVRGGRCAVLSCLPKEEQIGKCSTRGRKCCRRKK"},
    {"name": "ACTH", 
     "sequence": "SYSMEHFRWGKPVGKKRRPVKVYPNGAEDES"},
    {"name": "ASIP", 
     "sequence": "KKVVRPRTPLSAPCVATRNSCKPPAPACCDPCASCQCRFFRSACSCRVLSLNC"},
    #{"name": "Jn-IV",
     # "sequence": "GSICLEPKVVGPCTAYFRRFYFDSETGKCTPFIYGGCEGNSYVDEKLHACRAICRA"}
]

# 肽链长度和步移
peptide_lengths = [6, 7, 8, 9, 10]
peptide_offsets = [1, 1, 1, 1, 1]

# 输出文件名
output_file = r"C:\Users\123\Desktop\jobwork\subject\MC1R\screen\v1\text1.csv"

# 构建重叠肽库
overlap_library = []
for peptide in peptide_sequences:
    peptide_library = []  # 存储当前肽链的肽库
    for length in peptide_lengths:
        sequence_library = []  # 存储特定长度的肽库
        for offset in peptide_offsets:
            for i in range(0, len(peptide['sequence']) - length + 1, offset):
                subsequence = peptide['sequence'][i:i+length]
                sequence_library.append(subsequence)
        peptide_library.append((f"{peptide['name']} (Length {length})", sequence_library))
    overlap_library.extend(peptide_library)
    
pprint(overlap_library)

# 输出到CSV文件（覆盖旧文件）
with open(output_file, 'w+', newline='') as csvfile:
    writer = csv.writer(csvfile)

    # 按照不同肽链和长度输出到不同行和列
    for peptide_name, sequence_library in overlap_library:
        unique_sequences = list(OrderedDict.fromkeys(sequence_library))    
        #writer.writerow([peptide_name])
        writer.writerows([[peptide] for peptide in unique_sequences])
        #writer.writerow([])  # 添加空行以分隔不同序列的肽库
        
pprint(unique_sequences)