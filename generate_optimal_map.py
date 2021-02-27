from collections import defaultdict
import csv
from typing import Counter, Dict

import python_codon_tables


def read_csv_file(file_path):
    with open(file_path) as f:
        reader = csv.reader(f)
        records = list(reader)
    return records


def generate_codon_mapping_from_aa_mapping(codons_table) -> Dict:
    substitutions = {}
    for codon in codons_table.keys():
        frequency_table = codons_table[codon]
        codons = sorted(frequency_table)
        max_frequency = 0
        most_frequent_codon = None
        for c1 in codons:
            frequency = frequency_table[c1]
            if frequency > max_frequency:
                max_frequency = frequency
                most_frequent_codon = c1
        for c1 in codons:
            substitutions[codon] = most_frequent_codon
    return substitutions


vir_vac_combos = read_csv_file("side-by-side.csv")[1:]
# Generate a default dict with each amino acid and a list the most common swaps from the vaccine
most_frequent_codon_map = defaultdict(Counter)
for _, k, v in vir_vac_combos:
    most_frequent_codon_map[k][v] += 1


if __name__ == "__main__":
    substitutions = generate_codon_mapping_from_aa_mapping(most_frequent_codon_map)
    virvac = read_csv_file("side-by-side.csv")[1:]

    print(substitutions)

    matches = 0
    for element in virvac:
        _, vir, vac = element
        our = substitutions.get(
            vir, vir
        )  # stop codons are not encoded, just return stop codon
        if vac == our:
            matches += 1

    print("{:.1f}%".format(100 * matches / len(virvac)))
