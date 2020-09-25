#! /usr/bin/env python

import os
import subprocess
import edlib
import networkx as nx
import numpy as np
from . import utils
from .my_class import Monomers, Monomer 


def monomer_graph_analysis_check(input_fasta_file, monomer_fasta_file, output_prefix, 
    mean_monomer_len, head_to_tail_dist, min_fasta_len, thres_list):

    header_items = ("RID",
                    "Regularity",
                    "Read_Len",
                    "Thresh",
                    "#All_Monomers(clustered+not_clustered)",
                    "#mono_in_a_cluster",
                    "Isolates_(unclustered_monomers)",
                    "Clustered_monomer_fraction_in_read",
                    "#total_clusters_(distinct_monomers_in_HOR)",
                    "mean_identity_within_clusters",
                    "mean_identity_between_clusters",
                    "min_monomeric_period",
                    "max_monomeric_period",
                    "median_monomeric_period",
                    "Normalized_min_monomeric_period",
                    "Normalized_max_monomeric_period",
                    "min_head_to_tail_interval",
                    "max_head_to_tail_interval",
                    "median_head_to_tail_interval"
                    )
    header = "\t".join(header_items)
    monomer_db = {}

    seq_db = {}
    rc_map = dict(zip("ACGTacgtNn", "TGCAtgcaNn"))
    stats_formatting = "".join((
        "%d\t%f\t%d\t%d\t%d\t%f\t%d\t%f\t%f\t",
        "%d\t%d\t%d\t%f\t%f\t%d\t%d\t%d\n"))
    # output_prefix = os.path.basename(pread_filename).split(".")[0]
    regular_HORs_file = open(output_prefix + "_regularHORs.fa", 'w')
    inversions_file = open(output_prefix + "_inversions.fa", 'w')
    irregular_HORs_file = open(output_prefix + "_irregularHORs.fa", 'w')
    too_short_reads_file = open(output_prefix + "_too_short_reads.fa", 'w')
    no_HOR_reads_file = open(output_prefix + "_no_HOR_reads.fa", 'w')
    missing_monomer_file = open(output_prefix + "_missing_monomer.fa", 'w')
    regular_pattern_file = open(output_prefix + "_regularHORs_pattern.txt", 'w')
    irregular_pattern_file = \
        open(output_prefix + "_irregularHORs_pattern.txt", 'w')
    inversions_pattern_file = open(output_prefix + "_inversions_pattern.txt", 'w')
    stats_file = open(output_prefix + "_stats.txt", 'w')
    stats_file.write(header + "\n")

    # Print parameters
    print("Average monomer length: ", mean_monomer_len)
    print("Max head-to-tail distance: ", head_to_tail_dist)
    print("Shortest read length: ", min_fasta_len)
    print("Clustering thresolds: ", thres_list)


    monomers_db = {}
    seq_db = {}
    # IMPORT FASTA FILES #
    with open(input_fasta_file, 'r') as hin:
        for name, seq, qual in utils.readfq(hin):
            if len(seq) < min_fasta_len:
                print(">%s\n%s" % (name, seq), file = too_short_reads_file)
                continue
            # monomers_db[name] = Monomers(name, len(seq))
            seq_db[name] = seq

    
    # Load all monomers found in preads into monomer_db.
    # monomer_db[Read_ID] = [(start, end), sequence]
    with open(monomer_fasta_file, 'r') as hin:
        for name, seq, qual in utils.readfq(hin):
            # Parse the read tag.
            # Tag Format: ReadID/RangeStart_RangeEnd/Orientation
            rid, rng, orientation = name.split("/")
            # Skip if the read doesn't have any monomers.
            if rid not in seq_db:
                continue
            if rid not in monomers_db:
                monomers_db[rid] = Monomers(rid, len(seq_db[rid]))
            rng = rng.split("_")
            rng = int(rng[0]), int(rng[1])
            monomers_db[rid].add_monomer(Monomer(rng[0], rng[1], seq, orientation))
            # monomer_db.setdefault(rid, [])
            # monomer_db[rid].append((rng, seq, orientation))

    print(len(seq_db), " sequences read." , "Reads with monomers:", len(monomer_db.keys()))

    # RUN OVER ALL READS THAT CONTAIN MONOMERS #
    for rid, monomers in monomers_db.items():

        # import pdb; pdb.set_trace()
        monomers.sort_monomer()
        monomers.set_monomer_stat()
        monomers.comp_monomers()
        monomers.inversion_check()

        for threshold in thres_list:
            monomers.cluster_monomers(threshold)
            monomers.set_head_to_tail_stat()
            monomers.set_monomer_period_stat()

            """
            is_regular = False
            inversion_detected = False
            missing_monomer = False
            """

            # Exit the threshold loop if regularity, inversion or 
            # missing monomer detected.
            if inversion_detected:
                break
            elif ((max_abs_head_to_tail <= allowed_max_head_to_tail) and
                    (isolate_count == 0) and
                    (normalized_max_monomer_period <= 1.05) and
                    (normalized_min_monomer_period >= 0.95) and
                    (not inversion_detected)):
                # Mark as regular
                is_regular = True
                break
            elif isolate_count == 0 and \
              max_head_to_tail > 0.9 * avg_monomer_len:
                missing_monomer = True
                break           
        # END OF THE THRESHOLD LOOP #
        """
