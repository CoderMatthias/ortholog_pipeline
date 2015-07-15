#!/usr/bin/python3
'''Script reformats the output from synteny_checker.py and rearranges the output to rank by best predicted ortholog to worse predicted ortholog'''
import sys
import argparse

def argument_parser():
    '''Parses the given arguments from the input'''
    parser = argparse.ArgumentParser(description='Checks the detected genes to see if they are syntenic matches')
    parser.add_argument('input_', action='store', help='output (either ambiguous or unambiguous) from ortho_analysis.py')
    parser.add_argument('-l', action='store_true', dest='logs', help='write STDOUT and STDERR logs')
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

def combine_into_list (list_, header):
    """Opens an output file from synteny_checker.py, reformats, and adds to list which will be all combined synteny files"""
    formatted_list = [line[:header.index('syntenic')+1] + line[header.index('blastn_match'):] +\
                      line[header.index('syntenic')+1:header.index('blastn_match')] for line in list_]
    header = header[:header.index('syntenic')+1] + header[header.index('blastn_match'):] +\
             header[header.index('syntenic')+1:header.index('blastn_match')]
    return formatted_list, header

def sort_combined_list(list_, header):
    """
    Sorts orthologs according to the strength of the ortholog match.
    The most important to least important criteria are: synteny match, blastn, blastp
    """
    return sorted(list_, key=lambda x: (x[header.index('syntenic')] , x[header.index('blastn_match')] , x[header.index('blastp_match')]), reverse=True)

def remove_false_ortho(list_, header):
    header_terms = list(set([item for item in header if 'match' in item or 'syntenic' in item]))
    for line in list_:
        if line.count('no') == len(header_terms):
            line[header.index('ortholog')] = 'null'
    return list_

def write_output_file(header , species_name , combined_list):
    """Writes the output for the full combined, sorted list"""
    output_file_name = '9_{}_ortholog_analysis.tsv'.format(species_name)
    print 'Output saved as:', output_file_name
    with open(output_file_name, 'w+') as output_file:
        output_file.write('\t'.join(header) + '\n')
        for line in combined_list:
            output_line = '{}\n'.format('\t'.join(line))
            output_file.write(output_line)

def main():
    arg, species_name = argument_parser()
    LoL, header = source_file_LoL(arg.input_)
    formatted_list, header = combine_into_list(LoL, header)
    formatted_list = sort_combined_list(formatted_list, header)
    formatted_list = remove_false_ortho(formatted_list, header)
    write_output_file(header, species_name, formatted_list) 

if __name__ == "__main__":
    main()
