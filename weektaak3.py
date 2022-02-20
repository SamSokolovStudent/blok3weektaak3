from codondict import codons
import sys
import pandas as pd
from collections import Counter
import matplotlib.pyplot as plt
import numpy as np


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
            codon_bias_dict = codon_parser(codon_count_dict)
            codon_plotter(codon_bias_dict, element)


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
    final_list = ([])
    for keys, values in codon_values.items():
        codon_bias = {}
        for amino in codons.values():
            codon_bias[amino] = {}
        for keys_2, values_2 in values.items():
            amino_acid = codons[keys_2]
            codon_bias[amino_acid][keys_2] = values_2
        test = codon_bias
        final_list.append(test)
    return final_list


def codon_plotter(codon_bias_list, file_name):
    for element in codon_bias_list:
        df = pd.DataFrame(element)
        df.plot.bar()
        plt.title(file_name)
        plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=-3)
        f = plt.figure()
        f.set_figwidth(20)
        f.set_figheight(10)
        plt.show()


if __name__ == '__main__':
    main()
