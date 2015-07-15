#!/usr/bin/python3

import sys
import time

def genes_only_gff (sys_argv):
    with open(sys.argv[1], 'r') as source_file:
        name = source_file.name
        name = name[:3]
        print name
        gene_only_gff = []    
        line_list = source_file.read().split('\n')
        print 'gff read: DONE'
        for line in line_list:
            try:
                if line[0].strip() != '' and line[0] not in ('>' , '#' , 'A' , 'T' , 'C' , 'G' , 'N'):
                    data = line.strip().split('\t')
                    try:
                        if data[2] == 'gene':
                            data = [data[0]] + data[3:5] + data[8].split(';')
                            gene_only_gff.append(data)
                    except IndexError:
                        print 'INDEX ERROR ', data
                        continue 
            except IndexError:
                pass

        print 'Gene count = ', len(gene_only_gff)
        print 'genes_only_gff: DONE'
    return gene_only_gff , name

def write_output(gene_only_gff_list, name):
    with open('{}_gene_only.gff'.format(name), 'w+') as output_file:
        for gene in gene_only_gff_list:
            output_line = '\t'.join(gene)
            output_file.write(output_line + '\n')

def main():
    sys.stderr = open('ERROR.txt', 'w')
    sys.stdout = open('OUT.txt', 'w')
    st = time.time()
    gene_only_gff_list , name = genes_only_gff(sys.argv[1])
    write_output(gene_only_gff_list, name)

    print time.time() - st

if __name__ == '__main__':
    main()
