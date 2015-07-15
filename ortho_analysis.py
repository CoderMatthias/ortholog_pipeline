#!/usr/bin/python3

import sys
import argparse

def argument_parser():
    parser = argparse.ArgumentParser(description='Analyzed the combined blast results')
    parser.add_argument('input_', action='store', help='File of all the blast results lined up with D.mel ORs')
    parser.add_argument('-l', action='store_true', dest='logs', help='write STDOUT and STDERR logs')
    arg = parser.parse_args()
    spec_name = arg.input_.split('_')[1]
    return arg, spec_name

def logs(name, logs):
    if logs == True:
        sys.stdout = open('7_{}_analyzed_o.tsv'.format(name), 'w')
        sys.stderr = open('7_{}_analyzed_e.tsv'.format(name), 'w')

def source_file_LoL (input_, header='null'):
    '''takes file and turns it into a list of lists (LoL)'''
    with open(input_, 'r') as source_file:
        source_name = source_file.name
        print
        print 'Input file was:', source_name
        source_name = source_name[2:]
        line_list = [line for line in source_file.read().split('\n') if line.strip()]
        if line_list[0].startswith('#'):
            header = line_list[0].split('\t')
    print '\n{}'.format('-' * 25)
    print 'Total OR genes= {}\n'.format(len([line for line in line_list if not line.startswith('#')]))
    return [line.split('\t') for line in line_list if not line.startswith('#')], source_name, header

def insert_species_ortholog_column(list_, header):
    header.insert(header.index('FBgn_number') + 1, 'ortholog')
    for line in list_:
        line.insert(header.index('ortholog'), 'null')
    return list_, header

def parse_ambig_and_unambig(ortholog_list, header):
    unambig, ambig = [], []
    if 'blastn' in header:
        header.append('blastn_match')
        ortholog_list , unambig = one_matching_nucl_blast(ortholog_list, unambig, header.index('blastn'))
        ortholog_list , unambig , ambig = one_major_matching_nucl_blast(ortholog_list , unambig , ambig, header.index('blastn'))
    if 'blastp' in header:
        header.append('blastp_match')
        unambig , ambig = one_major_matching_prot_blast(ortholog_list + ambig, unambig, header.index('blastp'))
    return unambig, ambig, header

def one_matching_nucl_blast(ortholog_list , unambig_list, index):
    trim_ortho_list = []
    for line in ortholog_list:
        if line[index].strip() and not line[index + 1].strip() and float(line[index].split(',')[1]) > 200:
            line[line.index('null')] = line[index].split(',')[0]
            unambig_list.append(line)
        else:
            trim_ortho_list.append(line)
    print
    print 'Only one blastn match assessment'
    print 'unambiguous= {}'.format(len(unambig_list))
    print 'remain orthos= {}'.format(len(trim_ortho_list))
    return trim_ortho_list, unambig_list

def one_major_matching_nucl_blast(ortholog_list , unambig_list , ambig_list, index):
    trim_ortho_list = []
    for line in ortholog_list:
        try:
            top_hit, comparable = line[index].split(',')[1], line[index + 1].split(',')[1]
            if float(comparable) < float(top_hit) * 0.6 and float(top_hit) > 200:
                line[line.index('null')] = line[index].split(',')[0]
                unambig_list.append(line)
            else:
                ambig_list.append(line)
        except IndexError:
            trim_ortho_list.append(line)
    print
    print 'Only one major blastn match assessment'
    print 'unambiguous= {}'.format(len(unambig_list))
    print 'ambiguous= {}'.format(len(ambig_list))
    print 'remain orthos= {}'.format(len(trim_ortho_list))
    return trim_ortho_list, unambig_list, ambig_list

def one_major_matching_prot_blast(ortholog_list , unambig_list, index):
    ambig = []
    for line in ortholog_list:
        try:
            top_hit, comparable = line[index].split(',')[1], line[index + 1].split(',')[1]
            if float(comparable) < float(top_hit) * 0.6 and float(top_hit) > 200:
                line[line.index('null')] = line[index].split(',')[0]
                unambig_list.append(line)
            else:
                ambig.append(line)
        except:
            ambig.append(line)
    print
    print 'Only one major blastp match assessment'
    print 'unambiguous= {}'.format(len(unambig_list))
    print 'ambiguous= {}'.format(len(ambig))
    return unambig_list , ambig

def no_same_orthos(list_, header):
    dict_ = {line[0]: line for line in list_}
    list_ = sorted(list_, key=lambda x: x[header.index('ortholog')])
    assigned_orthos = [line[header.index('ortholog')] for line in list_]
    for i, line in enumerate(list_):
        if i == 0:
            continue
        if list_[i][header.index('ortholog')] == list_[i-1][header.index('ortholog')]:
            count, gene, same_orthos = -1, list_[i][header.index('ortholog')], []
            while list_[i + count][header.index('ortholog')] == gene:
                same_orthos.append(list_[i + count])
                count += 1
            same_orthos = sort_same_orthos_list(same_orthos, header)
            for line in same_orthos[1:]:
                possibles = [item for item in line[header.index('ortholog') + 1:] if item != '' and item != 'yes' and item != 'no']
                for item in possibles:
                    if item.split(',')[0] not in assigned_orthos:
                        line[header.index('ortholog')] = item.split(',')[0]
                        dict_[line[0]] = line
                        assigned_orthos.append(item.split(',')[0])
                        break
                else:
                    line[header.index('ortholog')] = 'null' 
                    dict_[line[0]] = line
    return [dict_[key] for key in dict_.keys()]

def sort_same_orthos_list(same_orthos, header):
    blast_list = list(sorted(set([item for item in header if 'blast' in item and 'match' not in item])))
    for i, item in enumerate(blast_list[:]):
        blast_list[i] = sorted([line for line in same_orthos if line[header.index(item)] != ''], key=lambda x: float(x[header.index(item)].split(',')[1]), reverse=True)
    blast_list = [item for sublist in blast_list for item in sublist]
    sorted_list, remove_redund = [], []
    for line in blast_list:
        if line not in remove_redund:
            sorted_list.append(line)
            remove_redund.append(line)
    return sorted_list

def separate_ambig_by_recip_blast(ambiguous_list, unambiguous_list, index, recip_index):
    remaining_orthos = []
    for line in ambiguous_list:
        if line[index].split(',')[0] == line[recip_index].split(',')[0]:
            unambiguous_list.append(line)
        else:
            remaining_orthos.append(line)
    print
    print 'Reciprical blast match assessment'
    print 'unambiguous= {}'.format(len(unambiguous_list))
    print 'ambiguous= {}'.format(len(remaining_orthos))
    return remaining_orthos, unambiguous_list

def blast_match_ortholog_check(list_, header):
    print '{}\n'.format('-' * 25)
    match_colm = [item[:item.index('_match')] for item in header if 'match' in item]
    for item in match_colm:
        for line in list_:
            if line[header.index(item)] != '':
                if line[header.index(item)].split(',')[0] == line[header.index('ortholog')]:
                    line.append('yes')
                else:
                    line.append('no')
            else:
                line.append('no')
    return list_

def fill_in_null_orthos(list_, header):
    blast_colm = list(sorted(set([item for item in header if 'blast' in item and 'match' not in item])))
    for line in list_:
        if line[header.index('ortholog')] == 'null':
            for item in blast_colm:
                if line[header.index(item)] != '':
                    line[header.index('ortholog')] = line[header.index(item)].split(',')[0]
                    break
    return list_ 

def write_output_file (output_file_name , list_to_write , header):
    '''writes the output file and prints name of the saved output file'''
    with open(output_file_name, 'w+') as output_file:
        print "Output saved as: {}".format(output_file.name)
        header = '\t'.join(header)
        output_file.write('{}\n'.format(header))
        for line in list_to_write:
            for item in line:
                if type(item) == list:
                    item = ','.join(map(str, item))
            output_line = '{}\n'.format('\t'.join(map(str, line)))
            output_file.write(output_line)

def main():
    arg, spec_name = argument_parser()
    logs(spec_name, arg.logs)
    OR_orthos, source_file_name, header = source_file_LoL(arg.input_)
    OR_orthos, header = insert_species_ortholog_column(OR_orthos, header)
    unambiguous, ambiguous, header = parse_ambig_and_unambig(OR_orthos, header)
    
    #ambiguous, unambiguous = separate_ambig_by_recip_blast(ambiguous, unambiguous, header.index('blastn'), header.index('rpblastn'))
    #ambiguous, unambiguous = separate_ambig_by_recip_blast(ambiguous, unambiguous, header.index('blastp'), header.index('rpblastp'))
    analyzed = unambiguous + ambiguous 
    analyzed = fill_in_null_orthos(analyzed, header)
    analyzed = no_same_orthos(analyzed, header)
    analyzed = blast_match_ortholog_check(analyzed, header)
    
    print
    print 'unambiguous where blastn top hit and blastp top hit don\'t match\tblastn\tblastp'
    count = 0
    for line in unambiguous:
        if line[header.index('blastn_match')] != line[header.index('blastp_match')]:
            if len(line[0]) > 6:
                print line[0], '\t', line[header.index('blastn')], '\t', line[header.index('blastp')]
            else:
                print line[0], '\t\t', line[header.index('blastn')], '\t', line[header.index('blastp')]
            count += 1
    print '# of unmatched', count
    print

    output_file_name = '7_{}_analyzed.tsv'.format(spec_name)
    write_output_file(output_file_name, analyzed, header)
    
    output_file_name = '7_{}_poss_dups.tsv'.format(spec_name)
    write_output_file(output_file_name, ambiguous, header)
    

if __name__ == '__main__':
    main()
