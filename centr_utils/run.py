#! /usr/bin/env python

from .split_to_monomer import split_to_monomer_check

def split_to_monomer_main(args):

    split_to_monomer_check(input_fasta_file, output_file, hmm_model_fwd, hmm_model_rev, min_monomer_len)


def monomer_graph_analysis_main(args):

    monomer_graph_analysis_check(input_fasta_file, monomer_fasta_file, output_prefix,
        mean_monomer_len, head_to_tail_dist, min_fasta_len, thres_list)


