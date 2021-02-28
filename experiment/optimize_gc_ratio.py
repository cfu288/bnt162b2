from collections import defaultdict
import csv
from typing import Dict

import python_codon_tables


def read_csv_file(file_path):
    with open(file_path) as f:
        reader = csv.reader(f)
        records = list(reader)
    return records


def replace_third_nucleotide_if_not_g_or_c(codon: str, reverse_map=dict) -> str:
    third_nucl = codon[2]
    if third_nucl == "G" or third_nucl == "C":
        return codon
    amino_acid = reverse_map.get(codon)
    new_codon_g = codon[:2] + "G"
    new_codon_c = codon[:2] + "C"

    if reverse_map.get(new_codon_c) == amino_acid:
        return new_codon_c
    if reverse_map.get(new_codon_g) == amino_acid:
        return new_codon_g
    return codon


def generate_codon_mapping_from_aa_mapping(species="h_sapiens_9606") -> Dict:
    substitutions = {}
    codons_table = python_codon_tables.get_codons_table(species)
    reverse_map = {}
    for aa, list_codons in codons_table.items():
        for codon in list_codons.keys():
            reverse_map[codon] = aa
    codons_to_map = reverse_map.keys()  # get all possible codos for a given aa
    for codon_to_map in codons_to_map:  # grab one of the codons
        best_value = replace_third_nucleotide_if_not_g_or_c(codon_to_map, reverse_map)
        substitutions[codon_to_map] = best_value

        if substitutions.get(codon_to_map) == None:  ##fallback if no match
            substitutions[codon_to_map] = codon_to_map
    return substitutions


if __name__ == "__main__":
    substitutions = generate_codon_mapping_from_aa_mapping()
    virvac = read_csv_file("side-by-side.csv")
    gc_ratio = 0.8
    total_nuc = 0
    gc_nuc = 0
    matches = 0
    with open("sample_solution_diff_1.txt", "w") as f:
        for element in virvac:
            _, vir, vac = element
            our = substitutions.get(
                vir, vir
            )  # stop codons are not encoded, just return stop codon
            if total_nuc != 0 and gc_nuc / total_nuc > gc_ratio:
                our = vir

            for char in our:
                if char == "G" or char == "C":
                    gc_nuc += 1
                total_nuc += 1

            if vac == our:
                matches += 1

        f.write(f"{_},{our},{vac}\n")
        print("{:.1f}%".format(100 * matches / len(virvac)))
