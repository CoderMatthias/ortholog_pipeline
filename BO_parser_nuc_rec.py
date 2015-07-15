#!/usr/bin/python3

import sys

def line_list_from_input (sys_argv):
    '''Open a file and make a line list with it's contents'''
    with open(sys_argv, 'r') as source_file:
        line_list = source_file.read().split('\n')
    source_file_name = source_file.name
    LoL = make_list_into_LoL(line_list)
    return LoL , source_file_name

def make_list_into_LoL (line_list):
    '''Take a line list and make a list of list (LoL)'''
    record , LoL = () , []
    for line in line_list:
        if not line.startswith('#') and line.strip() != '':
            record = line.split('\t')
            LoL.append(record)
    return LoL

def swap_mel_and_spec_columns (blast_LoL):
    '''Switches the dmel # and spec # columns, output in same format as non-reciprocal blast, which helps downstream in the pipeline'''
    swap_list = []
    for line in blast_LoL:
        line = [line[1]] + [line[0]] + line[2:]
        swap_list.append(line)
    return swap_list

def unique_list_by_column (list_to_filter , col_num):
    '''make a list of all the unique values in a given column'''
    unique_items = []
    for line in list_to_filter:
        if line[col_num] not in unique_items:
            unique_items.append(line[col_num])
    return unique_items

def make_blast_dict (blast_subset):
    '''make a dictionary of blast results where key = gene and value contains bitscore'''
    blast_dict = {}
    for gene in blast_subset:
        if gene[1] not in blast_dict:
            blast_dict[gene[1]] = [float(gene[-1])]
        else:
            blast_dict[gene[1]].append(float(gene[-1]))
    return blast_dict

def write_output(name , list_to_write):
    '''writes the output file and prints name of the saved output file'''
    with open('4_' + name[:-4] + '_parsed.tsv', 'w+') as output_file:
        print "Output saved as: {}".format(output_file.name)
        for line in list_to_write:
            output_line = '{}\n'.format('\t'.join(map(str, line)))
            output_file.write(output_line) 

def main():

    blastn_LoL , name = line_list_from_input(sys.argv[1])

    blastn_LoL = swap_mel_and_spec_columns(blastn_LoL) 

    mel_unique_genes = unique_list_by_column (blastn_LoL , 0)
    
    output_list = []
    for mel_gene in mel_unique_genes:
        blast_subset = []
        for gene in blastn_LoL:
            if gene[0] == mel_gene:
                blast_subset.append(gene)
    
        blast_dict = make_blast_dict (blast_subset)
        for gene in blast_subset:
            gene.append(sum(blast_dict[gene[1]]))
            if gene[:2] + [gene[-1]] not in output_list:
                output_list.append(gene[:2] + [gene[-1]])
    
    write_output(name, output_list)

if __name__ == '__main__':
    main()
