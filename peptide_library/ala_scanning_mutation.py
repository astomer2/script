import csv
from collections import OrderedDict
import openpyxl
from openpyxl import load_workbook



input_file = r"C:\Users\123\Desktop\jobwork\subject\MC1R\screen\v2\v2result.xlsx"
output_path = r"C:\Users\123\Desktop\jobwork\subject\MC1R\screen\v2\output.csv"

def ala_scan(seq):
  ala_seqs = []
  for i, aa in enumerate(seq):
    temp_seq = list(seq) 
    temp_seq[i] = 'A'
    ala_seqs.append(''.join(temp_seq))
  
  return ala_seqs


if input_file.endswith(".txt"):

  with open(input_file) as f:
    sequences = [line.strip() for line in f]

elif input_file.endswith(".xlsx"):

  wb = openpyxl.load_workbook(input_file)
  
  # 获取第二个工作表
  ws = wb.worksheets[1]  


  sequences = [cell.value for cell in ws['A']]
    

else:

  with open(input_file) as f:
    reader = csv.reader(f)
    sequences = [row[0] for row in reader]

with open(output_path, 'w') as out_f:
  writer = csv.writer(out_f)

  # 写入列头为原始序列
  writer.writerow([sequences]) 

  for seq in sequences:
    ala_seqs = ala_scan(seq)
    unique_sequences = list(OrderedDict.fromkeys(ala_seqs))

    # 每行写入原始序列
    row = [seq]

    # 对应Ala序列添加到行数据
    row.extend(unique_sequences) 

    writer.writerow(row)

print('完成')