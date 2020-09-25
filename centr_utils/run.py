#! /usr/bin/env python

from .split_to_monomer import split_to_monomer_check
from .monomer_graph_analysis import monomer_graph_analysis_check

def split_to_monomer_main(args):

    split_to_monomer_check(args.input_fasta_file, args.output_file, 
        args.hmm_model_fwd, args.hmm_model_rev, args.min_monomer_len)


def monomer_graph_analysis_main(args):

    monomer_graph_analysis_check(args.input_fasta_file, args.monomer_fasta_file, 
        args.output_prefix, args.mean_monomer_len, args.head_to_tail_dist, 
        args.min_fasta_len, args.thres_list) 


