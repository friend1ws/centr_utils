#! /usr/bin/env python

import os
import subprocess
from . import utils
from .monomer_class import Monomers, Monomer 


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
            seq_db[name] = seq

    
    # Load all monomers found in preads into monomer_db.
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

    print(len(seq_db), " sequences read." , "Reads with monomers:", len(monomer_db.keys()))

    # RUN OVER ALL READS THAT CONTAIN MONOMERS #
    for rid, monomers in monomers_db.items():

        print(rid)

        monomers.sort_monomer()
        monomers.set_monomer_stat()
        monomers.comp_monomers()

        for threshold in thres_list:
            monomers.cluster_monomers(threshold)
    
            if monomers.cluster_count > 1:
                monomers.set_head_to_tail_stat()
                monomers.set_monomer_period_stat()
                monomers.set_symbolic_pattern()

            monomers.state_check(head_to_tail_dist, mean_monomer_len)

            if monomers.inversion_detected == True:
                break
            elif monomers.is_regular == True:
                break
            elif monomers.missing_monomer == True:
                break
        # END OF THE THRESHOLD LOOP #

        if monomers.cluster_count > 1:
            fasta_tag = ">%s___%d__%d_%d__HOR%d" % \
                (rid, len(seq_db[rid]), monomers.HOR_start, monomers.HOR_end, int(monomers.cluster_count))
            fasta_seq = seq_db[rid][monomers.HOR_start:monomers.HOR_end + 1]
    
            stats_formatting = "%d\t%f\t%d\t%d\t%d\t%f\t%d\t%f\t%f\t%d\t%d\t%d\t%f\t%f\t%d\t%d\t%d"

            stats = stats_formatting % \
                (monomers.original_seq_len, threshold, len(monomers.monomer_list), 
                monomers.clustered_monomer_count, monomers.isolated_monomer_count, 
                float(monomers.clustered_monomer_len) / monomers.original_seq_len, 
                monomers.cluster_count, monomers.mean_identity_within_clusters, monomers.mean_identity_between_clusters,
                monomers.min_monomer_period, monomers.max_monomer_period, monomers.median_monomer_period,
                monomers.normalized_min_monomer_period, monomers.normalized_max_monomer_period, 
                monomers.min_head_to_tail, monomers.max_head_to_tail, monomers.median_head_to_tail)

        else:
            fasta_tag = ">" + rid
            fasta_seq = seq_db[rid]
            stats = "%d %f %d %e %d %s %e %s" % \
                (monomers.original_seq_len, threshold, len(monomers.monomer_list),
                monomers.clustered_monomer_count, len(monomers.monomer_list), ".", 
                monomers.cluster_count, "-1 "*10)

        if monomers.inversion_detected:
            print(rid + " V " + stats, file = stats_file)
            print(fasta_tag, file = inversions_file)
            print(fasta_seq, file = inversions_file)
            print(fasta_tag, file = inversions_pattern_file)
            print(monomers.symbolic_pattern, file = inversions_pattern_file)
        elif monomers.is_regular:
            print(rid + " R " + stats, file = stats_file)
            print(fasta_tag, file = regular_HORs_file)
            print(fasta_seq, file = regular_HORs_file)
            print(fasta_tag, file = regular_pattern_file)
            print(monomers.symbolic_pattern, file = regular_pattern_file)
        elif monomers.missing_monomer:
            print(rid + " M " + stats, file = stats_file)
            print(fasta_tag, file = missing_monomer_file)
            print(fasta_seq, file = missing_monomer_file)
        elif monomers.cluster_count > 1:
            print(rid + " I " + stats, file = stats_file)
            print(fasta_tag, file = irregular_HORs_file)
            print(fasta_seq, file = irregular_HORs_file)
            print(fasta_tag, file = irregular_pattern_file)
            print(monomers.symbolic_pattern, file = irregular_pattern_file)
        else:
            print(rid + " N " + stats, file = stats_file)
            print(fasta_tag, file = no_HOR_reads_file)
            print(fasta_seq, file = no_HOR_reads_file)


    regular_HORs_file.close()
    irregular_HORs_file.close()
    inversions_file.close()
    no_HOR_reads_file.close()
    too_short_reads_file.close()
    regular_HORs_file.close()
    missing_monomer_file.close()
    regular_pattern_file.close()
    irregular_pattern_file.close()
    inversions_pattern_file.close()
    stats_file.close()

