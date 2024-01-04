from dataclasses import dataclass
import pandas as pd
from pprint import pprint

@dataclass
class MutationGenerator:
    input_file: str
    output_file: str
    mutation_pos_list: list

def mutate_peptide(peptide, mutation_pos):
    """
    Generates a list of mutated peptides by replacing a specific amino acid in a given peptide.

    Parameters:
    - peptide (str): The original peptide string.
    - mutation_pos (int): The index position of the amino acid to be mutated.

    Returns:
    - mutated_peptides (List[str]): A list of mutated peptide strings.

    Example:
    >>> mutate_peptide("ABCDEF", 2)
    ['ABCAEF', 'ABCBEF', 'ABCCCF', 'ABCDFF', 'ABCEEF']
    """
    mutated_peptides = []
    for aa in "ACDEFGHIKLMNPQRSTVWY":
        if aa == peptide[mutation_pos]:
            continue
        mutated = peptide[:mutation_pos] + aa + peptide[mutation_pos+1:]
        mutated_peptides.append(mutated)
    return mutated_peptides

def read_peptide_seqs(file_path):
    """
    Read peptide sequences from a file.

    Parameters:
    file_path (str): The path to the file containing peptide sequences.

    Returns:
    list: A list of peptide sequences read from the file.

    Raises:
    None
    """
    if file_path.endswith(".csv"):
        df = pd.read_csv(file_path)
        return df["Sequence"].tolist()
    elif file_path.endswith(".txt"):
        with open(file_path) as f:
            return [line.strip() for line in f]
    else:
        print("Unsupported file format!")

def generate_mutated_library(peptide_seqs, mutation_pos_list):
    """
    Generates a mutated library of peptides based on the given peptide sequences
    and mutation positions.

    Parameters:
    - peptide_seqs (list): A list of peptide sequences.
    - mutation_pos_list (list): A list of mutation positions.

    Returns:
    - mutated_library (dict): A dictionary representing the mutated library of peptides,
      where the keys are sheet names and the values are the mutated peptides.

    """
    mutated_library = {}
    for i, peptide in enumerate(peptide_seqs):
        mutated_peptides = []
        for pos in mutation_pos_list:
            mutated = mutate_peptide(peptide, pos) 
            mutated_peptides.extend(mutated)
        mutated_library[f"Sheet{i+1}"] = mutated_peptides
    return mutated_library

def write_to_excel(output_file, mutated_library):
    """
        Writes the contents of a mutated library to an Excel file.

        Parameters:
            output_file (str): The path to the output Excel file.
            mutated_library (dict): A dictionary representing the mutated library, where the keys are peptide names and the values are lists of mutated peptides.

        Returns:
            None
    """
    writer = pd.ExcelWriter(output_file, engine='openpyxl') 
    df = pd.DataFrame()
    for i, (_, mutated_peptides) in enumerate(mutated_library.items()):
        col_name = f"Peptide {i+1}"
        df[col_name] = mutated_peptides
    df.to_excel(writer, index=False) 
    writer.close()

def main():
    input_file = r"D:\jobwork\subject\TRPV1\to_CaMK2\v1\ADCP\positive.txt"
    output_file = r'C:\Users\123\Downloads\mutated_library.xlsx'
    mutation_pos_list = [1, 2]

    mutation_generator = MutationGenerator(input_file, output_file, mutation_pos_list)

    # Testing
    peptide_seqs = read_peptide_seqs(mutation_generator.input_file) 
    mutated_library = generate_mutated_library(peptide_seqs, mutation_generator.mutation_pos_list)
    pprint(mutated_library)

    # Output to Excel
    write_to_excel(mutation_generator.output_file, mutated_library)
    print('完成')

if __name__ == "__main__":
    main()
