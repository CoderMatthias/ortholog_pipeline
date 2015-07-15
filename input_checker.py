#!/usr/bin/python2.7

import os
import argparse

def argument_parser():
    parser = argparse.ArgumentParser(description='Makes sure that all the necessary files are present to run the pipeline')
    parser.add_argument('files', nargs='+', action='store', help='output (either ambiguous or unambiguous) from ortho_analysis.py')
    return parser.parse_args()

def only_necessary_files(list_):
    begging_of_filenames = list(set([filename.split('_')[0] for filename in list_ if '_' in filename and 3 <= len(filename.split('_')[0]) <= 4 and not filename.endswith('.py')]))
    python_scripts = [filename for filename in list_ if filename.endswith('.py')]
    precomp_files = [filename for filename in list_ if filename.split('_')[0] in [item for item in begging_of_filenames if len(item) == 4]]
    spec_file_dict = {spec: [filename for filename in list_ if filename.startswith(spec)] for spec in begging_of_filenames if spec in [item for item in begging_of_filenames if len(item) ==3 and item != 'mel']}
    mel_files = [filename for filename in list_ if filename.startswith('mel')]
    return python_scripts, precomp_files, mel_files, spec_file_dict

def check_python_scripts(py_list):
    print
    needed_scripts = ('BO_parser_nuc.py', 'BO_parser_prot.py', 
                      'combine_synteny_output.py', 'make_dict.py', 
                      'ortho_analysis.py', 'synteny_checker.py', 
                      'BO_parser_nuc_rec.py', 'BO_parser_prot_rec.py', 
                      'gff_parser.py', 'make_mel_dict.py', 'ortho_list_w_blast_results.py')
    missing_scripts = [script for script in needed_scripts if script not in py_list]
    if len(missing_scripts) > 0:
        print 'You are missing or changed the name of the necessary script(s):'
        for script in missing_scripts:
            print script
    else:
        print 'All python scripts present'

def check_precomputed_files(precomp_list):
    print
    if 'gene_orthologs.tsv' not in precomp_list:
        ortho = 'gene_orthologs_fb_(release).tsv.gz'
    if 'fbgn' not in  [file_.split('_')[0] for file_ in precomp_list]:
        gn_tr_pp = 'fbgn_fbtr_fbpp_fb_(release).tsv.gz'
    if 'ortho' in locals() or 'gn_tr_pp' in locals():
        print 'You are missing the precomputed file(s):',
        if 'ortho' in locals():
            print ortho,
        if 'gn_tr_pp':
            print gn_tr_pp,
        print
        print 'You can get these from flybase under the releases section'
    else:
        print 'All precomputed files present'

def check_mel_files(mel_list):
    print
    if 'mel_gene_only.gff' not in mel_list or 'mel' not in [file_.split('_')[0] for file_ in mel_list if not file_.endswith('.gff')]:
        if 'mel.gff' in mel_list:
            print 'You need to run gff_parser.py on mel.gff'
        else:
            print 'You need the mel gff file'
            print 'You can download this from flybase under the genomes section'
        if 'mel' not in [file_.split('_')[0] for file_ in mel_list if not file_.endswith('.gff')]:
            print 'You need to name the list of the mel genes to check as something starting with mel'
    else:
        print 'All mel files present'

def check_drosophilia_species(spec_dict):
    print
    for spec in spec_dict.keys():
        needed_files = ('{}_all_prot.fasta'.format(spec),
                        '{}_blastn.tsv'.format(spec),
                        '{}_blastp.tsv'.format(spec),
                        '{}_gene_only.gff'.format(spec),
                        '{}_rpblastn.tsv'.format(spec),
                        '{}_rpblastp.tsv'.format(spec))
        for item in needed_files:
            if item not in spec_dict[spec]:
                nope = True
                if item == '{}_gene_only.gff'.format(spec) and '{}.gff' in spec_dict[spec]:
                    print 'You need to run gff_parser.py on {}.gff'.format(spec)
                else:
                    print 'You are missing file or need to change name to: {}'.format(item)
        if 'nope' not in locals():
            print 'All files are present for your species.'
            print 'Species you have files for =', spec_dict.keys()
        else:
            print 'You can find these files on flybase either under genomes section or the modencode section'
    print

def main():
    arg = argument_parser()
    python_scripts, precomp_files, mel_files, spec_file_dict = only_necessary_files(arg.files)
    check_python_scripts(python_scripts)
    check_precomputed_files(precomp_files)
    check_mel_files(mel_files)
    check_drosophilia_species(spec_file_dict)

if __name__ == '__main__':
    main()
