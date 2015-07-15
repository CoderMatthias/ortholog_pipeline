#!/usr/bin/python3

import sys
import argparse
import time

def argument_parser():
    parser = argparse.ArgumentParser(description='Format data to input into R')
    parser.add_argument('input_', nargs='+', action='store', help='Final analysis output file')
    parser.add_argument('-l', action='store_true', dest='logs', help='write STDOUT and STDERR logs')
    arg = parser.parse_args()
    return arg

def source_file_LoL(sys_argv):
    '''takes file and turns it into a list of lists (LoL)'''
    with open(sys_argv, 'r') as source_file:
        header = source_file.readline().strip().split('\t')
        species_name = source_file.name.split('_')[1]
        line_list = source_file.read().split('\n')
        LoL = [line.split('\t') for line in line_list if line.strip() and not line.startswith('#')]
    return LoL, header, species_name

def add_it_up(input_, syn_cong,):
    '''Meat of script, adds up various columns and gives combined numbers'''
    list_, header, species_name = source_file_LoL(input_)
    header_terms = [item for item in header if 'syntenic' in item or 'match' in item]
    compiled = []
    tot_ORs = len(list_)
    compiled.append(perfect_orthos(list_, species_name, header_terms))
    compiled.append(perfect_orthos_percentage(compiled , tot_ORs, species_name))
    compiled.append(perfect_orthos_percentage2(compiled , tot_ORs, species_name))
    for item in header_terms:
        compiled = add_yes_and_no(item, list_, tot_ORs, header.index(item), compiled, species_name, header)
    compiled.append(total_orthos(list_, species_name))
    compiled.append(total_orthos_percentage(compiled, tot_ORs, species_name))
    compiled.append(total_orthos_percentage2(compiled, tot_ORs, species_name))
    syn_cong = synteny_congruence(list_, syn_cong, species_name, header)
    return compiled, syn_cong

def perfect_orthos(LoL, species, header_terms):
    count = 0
    for line in LoL:
        if line.count('yes') == len(header_terms):
            count += 1
    return [species, count, 'perf ortho']

def perfect_orthos_percentage (list_ , total, name):
    return [name , float(list_[0][1])/total * 100 , 'perf ortho %']

def perfect_orthos_percentage2 (list_ , total, name):
    return [name , float(list_[0][1])/total * 100 , 'perf ortho %%']

def add_yes_and_no (header_term, LoL, tot_len, col_num, compiled, name, header):
    if '_' in header_term:
        col_num_2 = header.index(header_term[:header_term.index('_')])
    else:
        col_num_2 = 0
    ycount, value_present = 0, 0 
    for line in LoL:
        if line[col_num] == 'yes':
            ycount += 1 
        if len(line[col_num_2]) > 0:
            value_present += 1
    compiled.append([name , ycount , header_term])
    compiled.append(percentage(ycount, value_present, name, '{} %'.format(header_term)))
    compiled.append(percentage(ycount, tot_len, name, '{} %%'.format(header_term)))
    return compiled

def percentage (ycount , total, name, type_='null'):
    return [name , float(ycount)/total * 100 , type_]

def total_orthos (list_, species):
    count = 0 
    for line in list_:
        if 'yes' in line[3:10]:
            count += 1
    return [species , count , 'total ortho']

def total_orthos_percentage(list_, total, name):
    return [name , float(list_[-1][1])/total * 100 , 'total ortho %']

def total_orthos_percentage2 (list_ , total , name):
    return [name , float(list_[-2][1])/total * 100 , 'total ortho %%']

def synteny_congruence (list_, syn_cong, name, header, synteny_index=3):
    hi = [(header.index('blastn_match'), header.index('blastn')),
          (header.index('blastp_match'), header.index('blastp'))]
    for item in hi:
        y_count , n_count , b_count = 0 , 0 , 0
        for line in list_:
            if line[item[0]] == line[synteny_index] and line[item[1]].strip() != '':
                y_count += 1
            elif line[item[0]] != line[synteny_index] and line[item[1]].strip() != '':
                n_count += 1
            else:
                b_count += 1
        syn_cong.append([name, y_count, 'yes', '{}'.format(header[item[0]])])
        syn_cong.append([name, n_count, 'no', '{}'.format(header[item[0]])])
        syn_cong.append([name, b_count, 'null_blast', '{}'.format(header[item[0]])])
    return syn_cong

def syn_cong_percent (list_):
    return [[float(item)/202*100 if type(item) == int else item for item in line] for line in list_] 

def write_output_file (list_to_write , name='null'):
    output_file_name = 'R_{}_all_species_ortho_summary.tsv'.format(name)
    with open(output_file_name, 'w+') as output_file:
        for line in list_to_write:
            output_line = '{}\n'.format('\t'.join(map(str, line)))
            output_file.write(output_line)
    print 'Output saved as:', output_file_name

def main():
    syn_cong = []
    arg = argument_parser()
    combo = []
    for input_ in arg.input_:
        for_combo , syn_cong = add_it_up(input_, syn_cong)
        [combo.append(line) for line in for_combo]
    syn_cong.insert(0, ['species' , 'percent' , 'yes_no' , 'type'])
    syn_cong_perc = syn_cong_percent(syn_cong)
    write_output_file(syn_cong_perc, '3colsyn')
    head = ['species', 'total', 'type']
    norm_percent = [line for line in combo if line[2].count('%') == 1]
    norm_percent.insert(0, head)
    write_output_file(norm_percent, 'norm_percent')
    tru_percent = [line for line in combo if line[2].count('%') == 2]
    for line in tru_percent:
        line[2] = line[2][:-1]
    tru_percent.insert(0, head)
    write_output_file (tru_percent, 'tru_percent')

if __name__ == '__main__':
    main()
