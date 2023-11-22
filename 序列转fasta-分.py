import csv
import os

def seq_to_fasta(input_file, output_dir):

  if not os.path.exists(output_dir):
    os.makedirs(output_dir)  

  if input_file.endswith('.csv'):

    with open(input_file) as in_file:

      reader = csv.DictReader(in_file)

      for i, row in enumerate(reader):
        sequences = []
        sequences.append(row[columns_index])# 假设序列在第一列中

        output_file = os.path.join(output_dir, f"{i}.fasta")

        with open(output_file, 'w+') as out_file:

          out_file.write(f'>{i}\n')
          out_file.write(f'{sequence}\n')

  elif input_file.endswith('.txt'):

    with open(input_file) as in_file:

      for i, line in enumerate(in_file):

        sequence = line.strip()  

        output_file = os.path.join(output_dir, f"{i}.fasta")

        with open(output_file, 'w+') as out_file:

          out_file.write(f'>{i}\n') 
          out_file.write(f'{sequence}\n')

  print('转换完成!')

if __name__ == '__main__':

  input_file = r"C:\Users\123\Desktop\jobwork\subject\MC1R\screen\v1\adcp\seq-1.txt"
  output_dir = r'C:\Users\123\Desktop\jobwork\subject\MC1R\screen\v1\text'
  columns_index = 0
  seq_to_fasta(input_file, output_dir)