#! /usr/bin/env bash

python3 cycle.py alpha_curated.fa > alpha_curated_ext.fa   

centr_utils split_to_monomer alpha_curated_ext.fa alpha_curated_ext_monomer.fa ../example/alpha.hmm ../example/alpha.rc.hmm 

python3 change_label.py alpha_curated_ext_monomer.fa HOR.list.txt > alpha_curated_proc.fa

