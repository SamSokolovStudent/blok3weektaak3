from codondict import codons
import sys
import pandas as pd
from collections import Counter
import matplotlib.pyplot as plt


def main():
    file_opener()


def file_opener():
    file_name = sys.argv[1:]
    file_name = ["hiv1cdsnucleotide.txt", "hiv2cdsnucleotide.txt"]
    for element in file_name:
        with open(element, "r") as file:
            _seq_dict = file_reader(file)
            split_dict = codon_splitter(_seq_dict)
            codon_count_dict = codon_counter(split_dict)
            print(codon_count_dict)
            codon_bias_dict = codon_parser(codon_count_dict)
            print(codon_bias_dict)
            codon_plotter(codon_bias_dict)


def file_reader(file):
    _seq_dict = {}
    header = ""
    for line in file:
        if line.startswith(">"):
            header = line.split(" ")
            header = header[0]
            _seq_dict[header] = ""
        else:
            _seq_dict[header] += line.rstrip()
    return _seq_dict


def codon_splitter(_seq_dict):
    split_dict = {}
    for key in _seq_dict:
        sequence = _seq_dict[key]
        n = 3
        split_codons = [sequence[i:i+n] for i in range(0, len(sequence), n)]
        split_dict[key] = split_codons
    return split_dict


def codon_counter(_seq_dict):
    for key in _seq_dict:
        codon_buffer_dict = Counter(_seq_dict[key])
        codon_buffer_dict = dict(codon_buffer_dict)
        _seq_dict[key] = codon_buffer_dict
    return _seq_dict


def codon_parser(codon_values):
    codon_bias = {}
    for key in codon_values:
        buffer_list = list(codon_values[key])
        buffer_dict = dict(codon_values[key])
        for element in buffer_list:
            if element in codons:
                for amino in codons.values():
                    codon_bias[amino] = {}
                for codon_name, count in buffer_dict.items():
                    amino_acid = codons[codon_name]
                    codon_bias[amino_acid][codon_name] = count
    return codon_bias


def codon_plotter(codon_bias):
    super_codon_list = []
    super_values_list = []
    super_amino_list = []
    for aminos, _codons in codon_bias.items():
        codon_buffer_list = ([])
        value_buffer_list = ([])
        amino_buffer_list = ([])
        amino_buffer_list.append(aminos)

        print(aminos)
        for keys, values in _codons.items():
            codon_buffer_list.append(keys)
            value_buffer_list.append(values)
            super_codon_list.append(codon_buffer_list)
            super_values_list.append(value_buffer_list)
        super_amino_list.append(amino_buffer_list)
    print(super_codon_list)
    print(super_values_list)
    print(super_amino_list)


if __name__ == '__main__':
    main()
