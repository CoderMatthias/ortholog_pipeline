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

def remove_redundancy (list_w_redund):
    '''make a list removing multiple instances of item'''
    unique_list = []
    for line in list_w_redund:
        if line not in unique_list:
            unique_list.append(line)
    return unique_list

def write_output(name , list_to_write):
    '''writes the output file and prints name of the saved output file'''
    with open('1_' + name[:-4] + '_parsed.tsv', 'w+') as output_file:
        print "Output saved as: {}".format(output_file.name)
        for line in list_to_write:
            output_line = '{}\n'.format('\t'.join(map(str, line)))
            output_file.write(output_line) 

def main():

    blastn_LoL , name = line_list_from_input(sys.argv[1])

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
            if gene[0:2] + [gene[-1]] not in output_list:
                output_list.append(gene[0:2] + [gene[-1]])
    
    write_output(name, output_list)

if __name__ == '__main__':
    main()
