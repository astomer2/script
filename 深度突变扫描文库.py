import pandas as pd
from pprint import pprint

input_file = r"D:\jobwork\subject\TRPV1\to_CaMK2\v1\ADCP\positive.txt"
output_file = r'C:\Users\123\Downloads\mutated_library.xlsx'
mutation_pos_list = [1 ,2 ]


def mutate_peptide(peptide, mutation_pos):
    """对指定位点进行单位点突变"""
    mutated_peptides = []
    for aa in "ACDEFGHIKLMNPQRSTVWY":
        if aa == peptide[mutation_pos]:
            continue
        mutated = peptide[:mutation_pos] + aa + peptide[mutation_pos+1:]
        mutated_peptides.append(mutated)
    return mutated_peptides

def read_peptide_seqs(file_path):
    """读取不同格式的文件,获取多肽序列"""
    if file_path.endswith(".csv"):
        df = pd.read_csv(file_path)
        return df["Sequence"].tolist()
    elif file_path.endswith(".txt"):
        with open(file_path) as f:
            return [line.strip() for line in f]
    else:
        print("Unsupported file format!")
        
def generate_mutated_library(peptide_seqs, mutation_pos_list):
    mutated_library = {}
    for i, peptide in enumerate(peptide_seqs):
        mutated_peptides = []
        for pos in mutation_pos_list:
            mutated = mutate_peptide(peptide, pos) 
            mutated_peptides.extend(mutated)
        mutated_library[f"Sheet{i+1}"] = mutated_peptides
        
    return mutated_library

# 测试
peptide_seqs = read_peptide_seqs(input_file) 
mutated_library = generate_mutated_library(peptide_seqs, mutation_pos_list)


pprint(mutated_library)

# 输出到Excel
writer = pd.ExcelWriter(output_file, engine='openpyxl') 
df = pd.DataFrame()
for i, (_, mutated_peptides) in enumerate(mutated_library.items()):
    col_name = f"Peptide {i+1}"
    df[col_name] = mutated_peptides
    
#for sheet_name, mutated_peptides in mutated_library.items():
#    df = pd.DataFrame(mutated_peptides, columns=['Mutated Peptide'])
#    df.to_excel(writer, sheet_name=sheet_name, index=False)
#writer = pd.ExcelWriter('mutated_library.xlsx', engine='openpyxl')
df.to_excel(writer, index=False) 
writer.close()




