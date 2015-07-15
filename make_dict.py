#!/usr/bin/python

import sys
import argparse

def argument_parser():
    parser = argparse.ArgumentParser(description='This makes a species prot # to species gene # specifically for the genomes from modeencode')
    parser.add_argument('spec_all_prot', action='store', help='species_all_prot file')
    arg = parser.parse_args()
    return arg

def source_file_gene_list(input_arg):
    with open(sys.argv[1], 'r') as source_file:
        line_list = [line[1:] for line in source_file.read().strip().split('\n') if line.startswith('>')]
    return line_list

def make_dict(list_):
    dict_ = {}
    for line in list_:
        try:
            key = line[:line.index(' type=')]
        except ValueError:
            key = line
            dict_[key] = key
            continue
        term = line[line.index('parent=')+7:line.index('parent=')+18]
        dict_[key] = term
    return dict_

def write_output_file(dict_):
    output_name = '2_spec_dict.tsv'
    with open(output_name, 'w+') as output_file:
        print 'Output saved as:', output_name
        for key, term in dict_.iteritems():
            output_line = '{}\t{}\n'.format(key, term)
            output_file.write(output_line)

if __name__ == '__main__':
    arg = argument_parser()
    gene_name_list = source_file_gene_list(arg.spec_all_prot)
    gene_prot_dict = make_dict(gene_name_list)
    write_output_file(gene_prot_dict)
