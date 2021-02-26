from collections import defaultdict
import csv
from typing import Dict

import python_codon_tables


def read_csv_file(file_path):
    with open(file_path) as f:
        reader = csv.reader(f)
        records = list(reader)
    return records


def get_first_item_in_list_that_starts_with(z: str, y: list) -> str:
    return next(x for x in y if z[0] == x[0])


def generate_codon_mapping_from_aa_mapping(
    aa_mapping: Dict, species="h_sapiens_9606"
) -> Dict:
    substitutions = {}
    codons_table = python_codon_tables.get_codons_table(species)
    for AA, ideal_codons in aa_mapping.items():  # given our mappings of AA to codons
        codons_to_map = codons_table[AA].keys()  # get all possible codos for a given aa
        for codon_to_map in codons_to_map:  # grab one of the codons
            best_value = get_first_item_in_list_that_starts_with(
                codon_to_map[0], ideal_codons
            )
            substitutions[codon_to_map] = best_value

            if substitutions.get(codon_to_map) == None:  ##fallback if no match
                substitutions[codon_to_map] = ideal_codons[0]
    return substitutions


# Article: Codon optimality is tha major determinat of mRNA stability - Vladimir Presnyk http://dx.doi.org/10.1016/j.cell.2015.02.029
most_stable_aa_list = read_csv_file("codon_stabilization_coefficent.txt")[1:]
# Generate a default dict with each amino acid and a list of its most stable form
most_stable_aa_map = defaultdict(list)
for k, v in most_stable_aa_list:
    most_stable_aa_map[k].append(v)


if __name__ == "__main__":
    substitutions = generate_codon_mapping_from_aa_mapping(most_stable_aa_map)
    virvac = read_csv_file("side-by-side.csv")[1:]

    matches = 0
    for element in virvac:
        _, vir, vac = element
        our = substitutions.get(
            vir, vir
        )  # stop codons are not encoded, just return stop codon
        if vac == our:
            matches += 1
    print("{:.1f}%".format(100 * matches / len(virvac)))
