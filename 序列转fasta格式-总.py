import csv

def seq_to_fasta(input_file, output_file):
  
  # 检查输入输出文件类型
  if not isinstance(input_file, str):
    raise TypeError('Input file must be a string')
  
  if not isinstance(output_file, str):  
    raise TypeError('Output file must be a string')

  # 读取 CSV
  if input_file.endswith('.csv'):
    
    with open(input_file, 'r') as in_file:
      
      reader = csv.DictReader(in_file)
      
      with open(output_file, 'w+') as out_file:
      
        for i, row in enumerate(reader):
          sequences = []
          sequences.append(row[columns_index])# 假设序列在第一列中
          out_file.write(f'>{i}\n') 
          out_file.write(f'{sequence}\n')

  # 读取文本文件  
  elif input_file.endswith('.txt'):
    
    with open(input_file, 'r') as in_file:
    
      with open(output_file, 'w+') as out_file:
      
        for i, line in enumerate(in_file):
        
          sequence = line.strip()
          
          out_file.write(f'>seq{i}\n')
          out_file.write(f'{sequence}\n')
      
  print('转换完成!')
       
if __name__ == '__main__':

  input_file = r"C:\Users\123\Desktop\jobwork\subject\PRL\sequence\sequence.txt"
  output_file = r"C:\Users\123\Desktop\jobwork\subject\PRL\sequence\sequence.fasta"
  columns_index = 0
  seq_to_fasta(input_file, output_file)



