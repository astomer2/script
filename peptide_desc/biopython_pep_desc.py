import csv
from Bio.SeqUtils.ProtParam import ProteinAnalysis

def calculate_properties(sequence, pH):
    protein = ProteinAnalysis(sequence)
    molecular_weight = protein.molecular_weight()
    isoelectric_point = protein.isoelectric_point()
    instability_index = protein.instability_index()
    gravy = protein.gravy()
    charge = protein.charge_at_pH(pH)
    secondary_structure = protein.secondary_structure_fraction()  # [helix, turn, sheet]

    return molecular_weight, isoelectric_point, instability_index, gravy, charge, secondary_structure

def process_file(input_file, output_file, pH):
    sequences = []

    # 读取TXT文件中的多肽序列
    with open(input_file, 'r') as file:
        for line in file:
            sequence = line.strip()
            if sequence:
                sequences.append(sequence)

    # 计算每个多肽序列的属性并写入CSV表格
    with open(output_file, 'w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(['Sequence', 'Molecular Weight', 'Isoelectric Point',
                         'Instability Index', 'GRAVY', 'Charge', 'Helix Fraction',
                         'Turn Fraction', 'Sheet Fraction'])

        for sequence in sequences:
            molecular_weight, isoelectric_point, instability_index, gravy, charge, secondary_structure = calculate_properties(sequence, pH)
            writer.writerow([sequence, molecular_weight, isoelectric_point,
                             instability_index, gravy, charge, secondary_structure[0],
                             secondary_structure[1], secondary_structure[2]])

    print("计算完成！结果已保存到CSV文件中。")

# 替换为你的输入和输出文件路径以及pH值
input_file_path = r"C:\Users\123\Desktop\jobwork\subject\IL-8\seq.txt"
output_file_path = r"C:\Users\123\Desktop\jobwork\subject\IL-8\output1.csv"
pH_value = 7.4

process_file(input_file_path, output_file_path, pH_value)
