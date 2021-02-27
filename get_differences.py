from collections import defaultdict
import csv
from typing import Dict

import python_codon_tables


def read_csv_file(file_path):
    with open(file_path) as f:
        reader = csv.reader(f)
        records = list(reader)
    return records


if __name__ == "__main__":
    virvac = read_csv_file("side-by-side.csv")[1:]
    with open("differences.txt", "w") as f:
        for element in virvac:
            _, vir, vac = element

            if vir != vac:
                f.write(f"{_},{vir},{vac}\n")
