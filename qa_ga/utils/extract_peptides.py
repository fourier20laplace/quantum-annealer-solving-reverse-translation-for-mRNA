import random


def read_fasta_sequences(file_path):
    with open(file_path, 'r') as file:
        lines = file.readlines()

    sequences = []
    current_sequence = []

    for line in lines:
        line = line.strip()
        if line.startswith('>'):
            if current_sequence:
                sequences.append(''.join(current_sequence))
                current_sequence = []
        else:
            current_sequence.append(line)

    # Append the last sequence if file does not end with '>'
    if current_sequence:
        sequences.append(''.join(current_sequence))

    return sequences


def extract_peptides(sequence, peptide_length=20):
    peptides = []
    for i in range(0, len(sequence) - peptide_length + 1, peptide_length):
        peptides.append(sequence[i:i + peptide_length])
    return peptides


def get_peptides(file_path='./peptide_txt.txt', peptide_length=20, num_peptides=32, sel_all=False):
    sequences = read_fasta_sequences(file_path)
    all_peptides = []

    for sequence in sequences:
        peptides = extract_peptides(sequence, peptide_length)
        all_peptides.extend(peptides)
    if sel_all:
        return all_peptides
    else:
        random.seed(42)
        random_peptides = random.sample(all_peptides, num_peptides)
        return random_peptides
