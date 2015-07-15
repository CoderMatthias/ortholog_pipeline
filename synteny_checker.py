#!/usr/bin/python3

import sys
import argparse

def argument_parser():
    '''Parses the given arguments from the input'''
    parser = argparse.ArgumentParser(description='Checks the detected genes to see if they are syntenic matches')
    parser.add_argument('input_', action='store', help='output (either ambiguous or unambiguous) from ortho_analysis.py')
    parser.add_argument('mel_gff', action='store', help='gff of all melanogaster genes')
    parser.add_argument('species_gff', action='store', help='gff of all species genes that have a D.mel ortholog')
    parser.add_argument('-r', action='store', dest='range_', default=10, type=int, help='sets # of genes on either side for syntenic region')
    parser.add_argument('-t', action='store', dest='threshold', default=2, type=int, help='sets the amount of syntenic matches that must occur')
    parser.add_argument('-l', action='store_true', dest='logs', help='write STDOUT and STDERR logs')
    parser.add_argument('-synlog', action='store_true', dest='synlog', help='write a log showing the synteny of each gene')
    arg = parser.parse_args()
    spec_name = arg.input_.split('_')[1]
    logs(spec_name, arg.logs)
    return arg, spec_name

def logs(name, logs):
    '''Checks to see if a stdout and stderr should be written, and creates logs if so'''
    if logs:
        sys.stdout = open('8_{}_synteny_o.tsv'.format(name), 'w')
        sys.stderr = open('8_{}_synteny_e.tsv'.format(name), 'w')

def source_file_LoL (sys_argv):
    '''takes file and turns it into a list of lists (LoL)'''
    with open(sys_argv, 'r') as source_file:
        header = source_file.readline().strip().split('\t')
        line_list = [line for line in  source_file.read().split('\n') if line.strip() and not line.startswith('#')]
    return [line.split('\t') for line in line_list], header

def spec_to_mel_dict(list_, mel_index, spec_index):
    '''make a dictionary with key species gn# and value mel gn# from ortho_analysis.py output file'''
    return {line[spec_index]: line[mel_index] for line in list_ if line[spec_index] != 'null'}

def genes_with_orthos(spec_gff_input, ortho_dict, species_name):
    '''Checks to see if the gene on the .gff has a D.mel ortholog'''
    spec_gff_input, null = source_file_LoL(spec_gff_input)
    file_, path = open_ortho_file_or_mel_dict(spec_gff_input, species_name)
    spec_gff = []
    for line in spec_gff_input:
        if line[3][3:] in ortho_dict.keys():
            spec_gff.append(line[:3] + [line[3][3:]] + [ortho_dict[line[3][3:]]])
            continue 
        if path == 'dict':
            for item in line:
                if 'MELA:FB' in item:
                    spec_gff.append(line[:3] + [line[3][3:]] + [item[item.index('FB'):]])
                    break
        elif path == 'ortho':
            if line[3][3:] in file_.keys():
                spec_gff.append(line[:3] + [line[3][3:]] + [file_[line[3][3:]]])
    if path == 'dict':
        spec_gff = replace_FBpp_with_FBgn(spec_gff, file_)
    return spec_gff

def open_ortho_file_or_mel_dict(spec_gff_input, species_name):
    if spec_gff_input[0][3].startswith('ID=FBgn'):
        ortho_LoL, header = source_file_LoL('gene_orthologs.tsv')
        return {line[5]: line[0] for line in ortho_LoL if 'D{}\\'.format(species_name) in line[6]}, 'ortho'
    else:
        return source_file_LoL('2_mel_dict.tsv')[0], 'dict'

def replace_FBpp_with_FBgn(list_, mel_dict_LoL):
    '''Replaces the D.mel FBpp# with the FBgn# for each of the genes'''
    mel_dict = {line[0]: line[1] for line in mel_dict_LoL}
    print 'These FBpp\'s do not have a corrosponding FBgn:'
    for line in list_:
        try:
            line[4] = mel_dict[line[4]]
        except KeyError:
            if line[4].startswith('FBpp'):
                print line[4]
            pass
    print
    return list_

def spec_syn_dict(ortho_dict_keys, spec_gff, range_):
    '''make a dictionary with keys being species gn# and values being syntenic region using the species gff file'''
    return {spec_gff[i][3]: [line[4] for line in spec_gff[i-range_:i+range_+1]] for i, line in enumerate(spec_gff) if line[3] in ortho_dict_keys}

def mel_syn_dict(ortho_dict_values, spec_gff, range_):
    '''make a dictionary with keys being species gn# and values being syntenic region using the species gff file'''
    return {spec_gff[i][3][3:]: [line[3][3:] for line in spec_gff[i-range_:i+range_+1]] for i, line in enumerate(spec_gff) if line[3][3:] in ortho_dict_values}

def rev_ortho_dict(ortho_dict):
    return {value: key for key, value in ortho_dict.iteritems()}

def syntenic_checker (species_syntenic_regions, dmel_syntenic_regions, ortho_dict, genes_on_either_side, syntenic_match_threshold, species_name):
    syntenic_logger = []
    synteny_match_count , high_synteny_match_count , gene_count = 0 , 0 , 0
    not_syntenic , syntenic = ['not_syntenic'] , ['syntenic']
    sorted_keys = list(sorted(dmel_syntenic_regions.keys()))
    for key in sorted_keys:
        gene_count += 1
        counter = 0
        dmel_genes , dspec_genes = [] , []
        for mel_gene, spec_gene in zip(dmel_syntenic_regions[key], species_syntenic_regions[ortho_dict[key]]):
            if counter == genes_on_either_side:
                counter += 1
                dmel_genes.append('mel_gene')
                dspec_genes.append('{}_gene'.format(species_name))
                continue
            dmel_genes.append(mel_gene)
            dspec_genes.append(spec_gene)
            counter += 1
        synteny_match = 0
        for gene in dmel_genes:
            if gene in dspec_genes:
                synteny_match += 1
        syntenic_logger.append('{} ---  syntenic match in {} instances'.format(key , synteny_match))
        if synteny_match >= syntenic_match_threshold:
            syntenic.append(key)
            synteny_match_count += 1
            if synteny_match > genes_on_either_side:
                high_synteny_match_count += 1
            syntenic_logger.append('------------SYNTENIC   MATCH-------------')
        else:
            not_syntenic.append(key)
            syntenic_logger.append('********** NOT SYNTENIC MATCH ***********')
        syntenic_logger.append('   '.join(dmel_genes))
        syntenic_logger.append('   '.join(dspec_genes))
        syntenic_logger.append('')
    try:
        synteny_percentage = float(synteny_match_count) / float(gene_count) * 100
    except:
        synteny_percentage = 0
    try:
        high_synteny_percentage = float(high_synteny_match_count) / float(gene_count) * 100
    except:
        high_synteny_percentage = 0
    syntenic_logger.insert(0, '''
% of orthologs which would be predicted by the syntenic region: {}%
% of orthologs which half of the gene orthologs match within syntenic region: {}%

predicted orthologs by synteny analysis: {}
discounted orthologs by synteny analysis: {}

Synteny data for every ortholog in input list
(Lines are: Dmel FBgn---# syntenic matches in region; Syntenic match or not; Dmel gene region; Dana gene region) 
'''.format(synteny_percentage , high_synteny_percentage , len(syntenic)-1 , len(not_syntenic)-1))


    syntenic_logger.insert(0, '''
Synteny Checker Log

Number of genes on either side for syntenic region: {}
Number of genes that must match to count as syntenic region: {}'''.format(genes_on_either_side , syntenic_match_threshold))
    return syntenic , not_syntenic , syntenic_logger

def format_output(ortho_input, synteny_list, header, yn):
    if 'syntenic' not in header:
        header.insert(header.index('ortholog') + 1, 'syntenic')
    list_ = []
    for line in ortho_input:
        if line[header.index('FBgn_number')] in synteny_list:
            line.insert(header.index('syntenic'), yn)
            list_.append(line)
    return list_

def add_lines_from_input_not_yet_included(ortho_input, list_, header):
    present_genes = [line[header.index('FBgn_number')] for line in list_]
    re_check = []
    for line in ortho_input:
        if line[header.index('FBgn_number')] not in present_genes:
            re_check.append(line)
            line.insert(header.index('syntenic'), 'no')
            list_.append(line)
    return list_, re_check

def incorporating_recheck_data(synteny, list_, header):
    dict_ = {line[header.index('FBgn_number')]: line for line in list_}
    print list_
    temp = []
    for line in synteny:
        if line[header.index('FBgn_number')] in dict_.keys():
            print line
            line = dict_[line[header.index('FBgn_number')]]
            print line
            temp.append(line)
        else:
            temp.append(line)
    return temp

def write_synteny_output (synteny_file, header, species_name):
    '''Write an output file of given name'''
    output_name = '8_{}_synteny_analyzed.tsv'.format(species_name) 
    print 'Output saved as:', output_name
    with open(output_name, 'w+') as output_file:
        output_file.write('{}\n'.format('\t'.join(header)))
        for line in synteny_file:
            output_line = '{}\n'.format('\t'.join(line))
            output_file.write(output_line)

def write_log_file(logger_list , species_name, synlog):
    '''Write the log file for the cycle iteration'''
    if synlog:
        logger_output_name = '8_{}_synteny_log.txt'.format(species_name)
        with open(logger_output_name, 'w+') as logger:
            logger.write('\n'.join(logger_list))

def main ():
    arg, species_name = argument_parser()
    ortho_input, header = source_file_LoL(arg.input_)
    mel_gff, null = source_file_LoL(arg.mel_gff)
    ortho_dict = spec_to_mel_dict(ortho_input, header.index('FBgn_number'), header.index('ortholog'))
    species_gff = genes_with_orthos(arg.species_gff, ortho_dict, species_name)
    spec_syn = spec_syn_dict(ortho_dict.keys(), species_gff, arg.range_)
    mel_syn = mel_syn_dict(ortho_dict.values(), mel_gff, arg.range_)
    rortho_dict = rev_ortho_dict(ortho_dict)
    syntenic, not_syntenic, syntenic_logger = syntenic_checker (spec_syn, mel_syn, rortho_dict, arg.range_, arg.threshold, species_name)
    synteny = format_output(ortho_input, syntenic, header, 'yes') + format_output(ortho_input, not_syntenic, header, 'no')
    synteny, recheck = add_lines_from_input_not_yet_included(ortho_input, synteny, header)
    write_log_file(syntenic_logger, species_name, arg.synlog)
    write_synteny_output(synteny, header, species_name)

if __name__ == '__main__':
    main()
