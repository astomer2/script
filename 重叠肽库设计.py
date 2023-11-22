import pandas as pd

# 肽链名称和序列
peptide_sequences = {"PRL": ["PICPGGAARCQVTLRDLFDRAVVLSHYIHNLSSEMFSEFDKRYTHGRGFITKAINSCHTSSLATPEDKEQAQQMNQKDFLSLIVSILRSWNEPLYHLVTEVRGMQEAPEAILSKAVEIEEQTKRLLEGMELIVSQVHPETKENEIYPVWSGLPSLQMADEESRLSAYYNLLHCLRRDSHKIDNYLKLLKCRIIHNNNC"],
                     "PRL-surface": ["PAQVDLRAVLHYHNITINPEDKEQEQREGLHRRHKDNYKLCNNC"]}

# 肽链长度和步移
peptide_lengths = [6, 7, 8, 9, 10]
peptide_offsets = [1, 1, 1, 1, 1]

# 创建空的数据框
df = pd.DataFrame(columns=["Peptide Name", "Peptide Sequence"])

# 生成肽链库
data_to_append = []  # 用于累积要添加的数据

for name, sequence in peptide_sequences.items():
    for length, offset in zip(peptide_lengths, peptide_offsets):
        for i in range(0, len(sequence[0]) - length , offset):
            peptide_name = f"{name}_{length}_{i}"
            peptide_seq = sequence[0][i:i + length]
            data_to_append.append({"Peptide Name": peptide_name, "Peptide Sequence": peptide_seq})

df = pd.concat([df, pd.DataFrame(data_to_append)], ignore_index=True)

# 保存到CSV文件
output_file = r"C:\Users\123\Desktop\jobwork\subject\PRL\overlap_library.csv"
df.to_csv(output_file, index=False)

print(f"Library has been saved to {output_file}")
