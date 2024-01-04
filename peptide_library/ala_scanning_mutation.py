import csv
from collections import OrderedDict
import openpyxl

def ala_scan(seq):
    """
    Generate alanine (Ala) sequences by replacing each amino acid in the given sequence with Ala.
    
    Parameters:
    - seq (str): Input amino acid sequence.

    Returns:
    - list: List of generated Ala sequences.
    """
    ala_seqs = []
    for i, aa in enumerate(seq):
        temp_seq = list(seq)
        temp_seq[i] = 'A'
        ala_seqs.append(''.join(temp_seq))
  
    return ala_seqs

def read_sequences_from_txt(file_path):
    """
    Read sequences from a text file.

    Parameters:
    - file_path (str): Path to the input text file.

    Returns:
    - list: List of sequences read from the file.
    """
    with open(file_path) as f:
        sequences = [line.strip() for line in f]
    return sequences

def read_sequences_from_excel(file_path):
    """
    Read sequences from an Excel file.

    Parameters:
    - file_path (str): Path to the input Excel file.

    Returns:
    - list: List of sequences read from the file.
    """
    wb = openpyxl.load_workbook(file_path)
    ws = wb.worksheets[1]  # Assuming you want the second worksheet
    sequences = [cell.value for cell in ws['A']]
    return sequences

def read_sequences_from_csv(file_path):
    """
    Read sequences from a CSV file.

    Parameters:
    - file_path (str): Path to the input CSV file.

    Returns:
    - list: List of sequences read from the file.
    """
    with open(file_path) as f:
        reader = csv.reader(f)
        sequences = [row[0] for row in reader]
    return sequences

def write_sequences_to_csv(output_path, header, data):
    """
    Write sequences to a CSV file.

    Parameters:
    - output_path (str): Path to the output CSV file.
    - header (list): List of column headers.
    - data (list): List of rows containing sequences.

    Returns:
    - None
    """
    with open(output_path, 'w', newline='') as out_f:
        writer = csv.writer(out_f)
        writer.writerow(header)
        writer.writerows(data)

def main(input_file, output_path):
    if input_file.endswith(".txt"):
        sequences = read_sequences_from_txt(input_file)
    elif input_file.endswith(".xlsx"):
        sequences = read_sequences_from_excel(input_file)
    else:
        sequences = read_sequences_from_csv(input_file)

    # Prepare data for writing to CSV
    header = ["Original Sequence", "Ala Sequences"]
    data = []

    for seq in sequences:
        ala_seqs = ala_scan(seq)
        unique_sequences = list(OrderedDict.fromkeys(ala_seqs))
        row = [seq, ', '.join(unique_sequences)]  # Combine Ala sequences into a comma-separated string
        data.append(row)

    write_sequences_to_csv(output_path, header, data)

if __name__ == "__main__":
    input_file = r"C:\Users\123\Desktop\jobwork\subject\MC1R\screen\v2\v2result.xlsx"
    output_path = r"C:\Users\123\Desktop\jobwork\subject\MC1R\screen\v2\output.csv"
    main(input_file, output_path)
    print('完成')