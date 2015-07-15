#!/bin/bash

echo -n "What is the species 3-letter abreviated Drosophila species name: "
read -e NAME

python BO_parser_nuc.py ${NAME}_blastn.tsv
python make_mel_dict.py fbgn_fbtr_fbpp_fb_2015_01.tsv
python make_dict.py ${NAME}_all_prot.fasta
python BO_parser_prot.py ${NAME}_blastp.tsv 2_mel_dict.tsv 2_spec_dict.tsv
python BO_parser_nuc_rec.py ${NAME}_rpblastn.tsv
python BO_parser_prot_rec.py ${NAME}_rpblastp.tsv 2_mel_dict.tsv 2_spec_dict.tsv
python ortho_list_w_blast_results.py mel_OR_genes.tsv -bn 1_${NAME}_blastn_parsed.tsv -brn 4_${NAME}_rpblastn_parsed.tsv -bp 3_${NAME}_blastp_parsed.tsv -brp 5_${NAME}_rpblastp_parsed.tsv
python ortho_analysis.py 6_${NAME}_blast_analyzed.tsv
python synteny_checker.py 7_${NAME}_analyzed.tsv mel_gene_only.gff ${NAME}_gene_only.gff -synlog
python combine_synteny_output.py 8_${NAME}_synteny_analyzed.tsv
