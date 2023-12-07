import pandas as pd

def read_fasta(file_path):
    """
    Read peptide sequences from a FASTA file.

    Parameters:
    - file_path (str): Path to the input FASTA file.

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
                    peptide_name = f"{name}_{length}_{i}"
                    peptide_seq = sequence[i:i + length]
                    data_to_append.append({"Peptide Name": peptide_name, "Peptide Sequence": peptide_seq})

    return pd.DataFrame(data_to_append)

def save_library_to_csv(data_frame, output_file):
    """
    Save DataFrame to a CSV file.

    Parameters:
    - data_frame (pd.DataFrame): DataFrame to be saved.
    - output_file (str): Path to the output CSV file.
    """
    data_frame.to_csv(output_file, index=False)
    print(f"Library has been saved to {output_file}")

def main():
    peptide_sequences = read_fasta(fasta_file)
    # Generate and save peptide library to CSV
    peptide_library_df = generate_peptide_library(peptide_sequences, peptide_lengths, peptide_offsets)
    save_library_to_csv(peptide_library_df, output_file)

if __name__ == "__main__":
    # Path to the FASTA file containing peptide sequences
    fasta_file = r"path/to/your/peptide_sequences.fasta"
    output_file = r"C:\Users\123\Desktop\jobwork\subject\PRL\overlap_library.csv"
    # Define peptide lengths and offsets
    peptide_lengths = [6, 7, 8, 9, 10]
    peptide_offsets = [1, 1, 1, 1, 1]