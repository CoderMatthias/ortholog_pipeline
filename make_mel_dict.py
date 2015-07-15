#!/usr/bin/python

import sys

source_file = open(sys.argv[1], 'r')
output_name = '2_mel_dict.tsv'
output_file = open(output_name, 'w+')

line_list = source_file.read().split('\n')

record , filer = () , []
for line in line_list:
    if not line.startswith('#') and line.strip() != '':
        record = line.split('\t')
        filer.append(record)

trim_ll = []
for line in filer:
    if len(line) == 3:
        trim_ll.append(line)

mel_dict = {}
for line in trim_ll:
    mel_dict[line[2]] = line[0]

for key, term in mel_dict.items():
    output_line = '%s\t%s\n' % (key , term)
    output_file.write(output_line)


print 'Output saved as:', output_name
source_file.close()
output_file.close()

