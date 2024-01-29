from collections import OrderedDict

def read_fasta(file_path):
    """
    Read peptide sequences from a FASTA file.

    Parameters:
    - file_path (str): Path to the input FASTA or TXT file.

    Returns:
    - dict: Dictionary with peptide names as keys and sequences as values.
    """
    peptide_sequences = {}
    with open(file_path, 'r') as f:
        lines = f.readlines()
        current_name = None
        current_sequence = []
        for line in lines:
            line = line.strip()
            if line.startswith(">"):
                if current_name is not None:
                    peptide_sequences[current_name] = [''.join(current_sequence)]
                current_name = line[1:]
                current_sequence = []
            else:
                current_sequence.append(line)
        if current_name is not None:
            peptide_sequences[current_name] = [''.join(current_sequence)]
    return peptide_sequences

def generate_peptide_library(peptide_sequences, lengths, offsets):
    """
    Generate a peptide library based on input sequences, lengths, and offsets.

    Parameters:
    - peptide_sequences (dict): Dictionary with peptide names as keys and sequences as values.
    - lengths (list): List of peptide lengths.
    - offsets (list): List of offsets.

    Returns:
    - pd.DataFrame: DataFrame containing the generated peptide library.
    """
    data_to_append = []
    for name, sequence_list in peptide_sequences.items():
        for length, offset in zip(lengths, offsets):
            for sequence in sequence_list:
                for i in range(0, len(sequence) - length + 1, offset):
                    peptide_seq = sequence[i:i + length]
                    data_to_append.append( peptide_seq)

    unique_sequences = list(OrderedDict.fromkeys(data_to_append))  

    return unique_sequences

def save_library_to_txt(unique_sequences, output_file):
    """
    Save DataFrame to a CSV file.

    Parameters:
    - data_frame (pd.DataFrame): DataFrame to be saved.
    - output_file (str): Path to the output CSV file.
    """
    with open(output_file, 'w', newline='') as f:
        for seq in unique_sequences:
            f.write(f"{seq}\n")


def main(fasta_file, output_file, peptide_lengths, peptide_offsets):
    peptide_sequences = read_fasta(fasta_file)
    # Generate and save peptide library to CSV
    unique_sequences = generate_peptide_library(peptide_sequences, peptide_lengths, peptide_offsets)
    save_library_to_txt(unique_sequences, output_file)

if __name__ == "__main__":
    # Path to the FASTA file containing peptide sequences
    fasta_file = '/mnt/nas1/lanwei-125/MC5R/Sequence/target_Sequence.txt'
    # 传入的fasta中可以是一条序列，也可以是多条，但是需要规范写法

    output_file = '/mnt/nas1/lanwei-125/MC5R/Sequence/split_sequence.txt'
    # Define peptide lengths and offsets
    peptide_lengths = [4, 5]
    peptide_offsets = [1, 1]
    main(fasta_file, output_file, peptide_lengths, peptide_offsets)