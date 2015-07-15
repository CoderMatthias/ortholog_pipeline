#!/usr/bin/python3

import sys
import time
'''
This python script takes the blastp tab-deliminated data and parses it to include unique hits and an additive bit score
Requires: output from blastp, a dictionary to convert FBpp#s to FBgn#s for both Dmen and other D species
Usage:   python BO_parser_prot.py blastp_output.tsv Dmel_FBpp_to_FBgn.tsv species_FBpp_to_FBgn.tsv
'''
start_time = time.time()

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

def FBpp_to_FBgn_dict (dict_sys_argv):
    '''Converts the FBpp <-> FBgn file and to a dictionary'''
    LoL , nullname = line_list_from_input (dict_sys_argv)

    out_dict = {}
    for line in LoL:
        out_dict[line[0]] = line[1]
    return out_dict

def swap_mel_and_spec_columns (blast_LoL):
    '''Switches the dmel # and spec # columns, output in same format as non-reciprocal blast, which helps downstream in the pipeline'''
    swap_list = []
    for line in blast_LoL:
        line = [line[1]] + [line[0]] + line[2:]
        swap_list.append(line)
    return swap_list

def replace_FBpp_w_FBgn (pp_to_gn_dict , list_to_switch , column_to_switch):
    '''Replaces protein number (FBpp) with the gene number (FBgn)'''
    for line in list_to_switch:
        line[column_to_switch] = pp_to_gn_dict[line[column_to_switch]]
    return list_to_switch

def column_value_unique_list (line_list , column_number):
    '''Make list of all unique items in column of list'''
    unique_list = []
    for line in line_list:
        if line[column_number] not in unique_list:
            unique_list.append(line[column_number])
    return unique_list

def make_blast_dict (blast_subset):
    '''make a dictionary of blast results where key = gene and value contains bitscore'''
    blast_dict = {}
    for line in blast_subset:
        if line[1] not in blast_dict:
            blast_dict[line[1]] = [float(line[-1])]
        else:
            blast_dict[line[1]].append(float(line[-1]))
    return blast_dict

def write_output(name , list_to_write):
    '''Write an output file from a list'''
    output_file_name = '5_{}_parsed.tsv'.format(name[:-4])
    print 'Output file saved as: {}'.format(output_file_name)
    with open(output_file_name, 'w') as output_file:
        for line in list_to_write:
            output_line = '\t'.join(map(str, line))
            output_file.write(output_line + '\n')

def main():
    blast_list , source_file_name = line_list_from_input (sys.argv[1])
    m_dict = FBpp_to_FBgn_dict (sys.argv[2])
    s_dict = FBpp_to_FBgn_dict (sys.argv[3])
    blast_list = swap_mel_and_spec_columns (blast_list)
    blast_list = replace_FBpp_w_FBgn (m_dict , blast_list , 0)
    blast_list = replace_FBpp_w_FBgn (s_dict , blast_list , 1)
    unique_mel_genes = column_value_unique_list (blast_list , 0)

    output_list = []
    for mel_gene in unique_mel_genes:
        blast_subset , new_blast_list = [] , []
        for line in blast_list:
            if line[0] == mel_gene:
                blast_subset.append(line)
            else:
                new_blast_list.append(line)
        blast_list = new_blast_list 
    
        blast_dict = make_blast_dict (blast_subset)
       
        for blast in blast_subset:
            blast.append(sum(blast_dict[blast[1]]))
            if blast[:2] + [blast[-1]] not in output_list:
                output_list.append(blast[:2] + [blast[-1]])
 
    write_output(source_file_name , output_list)
#    print (time.time()-start_time)

if __name__ == '__main__':
    main()    
