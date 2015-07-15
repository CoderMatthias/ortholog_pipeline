#!/usr/bin/python3

import sys
import argparse
import time

def argument_parser():
    parser = argparse.ArgumentParser(description='Combine the melanogaster blast results')
    parser.add_argument('input_', action='store', help='File of D. mel OR genes')
    parser.add_argument('-bn', action='store', dest='nuc', default='null', help='blastn output from BO_parser_nucl.py')
    parser.add_argument('-brn', action='store', dest='rpnuc', default='null', help='reciprical blastn output from BO_parser_nucl.py')
    parser.add_argument('-bx', action='store', dest='nucex', default='null', help='blastn_ex output from BO_parser_nucl.py')
    parser.add_argument('-bp', action='store', dest='prot', default='null', help='blastp output from BO_parser_prot.py')
    parser.add_argument('-brp', action='store', dest='rpprot', default='null', help='reciprical blastp output from BO_parser_prot.py')
    parser.add_argument('-nh', action='store', dest='numhits', default=5, type=int, help='top # of hits included in output file from each blast input')
    parser.add_argument('-l', action='store_true', dest='logs', help='write STDOUT and STDERR logs')
    arg = parser.parse_args()
    species_name =  arg.nuc.split('_')[1]
    return arg, species_name

def logs(name, logs):
    if logs == True:
        sys.stdout = open('6_{}_blast_analyzed_o.tsv'.format(name), 'w')
        sys.stderr = open('6_{}_blast_analyzed_e.tsv'.format(name), 'w')

def source_file_LoL (file_input):
    '''takes file and turns it into a list of lists (LoL)'''
    with open(file_input, 'r') as source_file:
        source_file_name = source_file.name
        line_list = [line for line in source_file.read().split('\n') if line.strip()]
        if line_list[0].startswith('#'):
            header = line_list[0].split('\t')
        else:
            header = 'null'
    return [line.split('\t') for line in line_list if not line.startswith('#')], source_file_name, header

def try_blast_files(blastn, rpblastn, blastnx, blastp, rpblastp, numhits, list_, header):
    blasts = [blastn, rpblastn, blastnx, blastp, rpblastp]
    for blast in blasts:
        if blast != 'null':
            list_, header = open_blast(list_, blast, numhits, header) 
    return list_, header

def open_blast (formatted_gene_list, sys_argv, blast_hits, header):
    '''opens the blast parsed data and combines it with formatted data'''

    blast_LoL , blast_name, null = source_file_LoL(sys_argv)
    blast_LoL = [[item[:item.index('-')] if '-' in item else item for item in line] for line in blast_LoL]
    species_name, blast_name = blast_name.split('_')[1], blast_name.split('_')[2]
    print
    print 'blast input', blast_name
    header += [blast_name] * blast_hits
    blast_dict = {key[0]: sorted([[line[1]] + [float(line[2])] for line in blast_LoL if line[header.index('#genes')] == key[0]], key = lambda x : x[1], reverse=True) for key in blast_LoL}
    blast_dict = {key: value[:blast_hits] if len(value) >= blast_hits else value + [''] * (blast_hits - len(value)) for key, value in blast_dict.iteritems()}
    formatted_gene_list = [line + blast_dict[line[header.index('FBgn_number')]] if line[header.index('FBgn_number')] in blast_dict.keys() else line + [''] * blast_hits for line in formatted_gene_list]
    print 'Completed open blast for:', blast_name
    return formatted_gene_list, header

def write_output_file(list_to_write, name, header):
    with open('6_' + name + '_blast_analyzed.tsv', 'w+') as output_file:
        print 'Output saved as:', output_file.name
        print
        output_file.write('{}\n'.format('\t'.join(header)))
        for line in list_to_write:
            output_line = [','.join(map(str, item)) if isinstance(item, list) else item for item in line]
            output_line = "{}\n".format('\t'.join(output_line))
            output_file.write(output_line)

def main():
    st = time.time()
    arg, species_name = argument_parser()
    logs(species_name, arg.logs)
    mel_ORs, null, header = source_file_LoL(arg.input_)
    mel_ORs, header = try_blast_files (arg.nuc, arg.rpnuc, arg.nucex, arg.prot, arg.rpprot, arg.numhits, mel_ORs, header)
    write_output_file(mel_ORs, species_name, header)

if __name__ == '__main__':
    main()
