#! /usr/bin/env python

import os
import subprocess
import edlib
import networkx as nx
import numpy as np
from . import utils
from class import Monomers, 
class Monomers(object):
    
    def __init__(self, origina_seq_len, original_seq_id):
        self.original_seq_len = None
        self.original_seq_id = NOne
        self.aln_info = []
        self.total_monomer_len = None
        self.is_regular = False
        self.inversion_detected = False
        self.missing_monomer = False
        self.monomer_list = []


class Monomer(object):
    
    def __init__(self):
        self.start = None
        self.end = None
        self.seq = None
        self.orientation = None
        self.cluster_id = None

class Cluster(ojbect):

    def __init__(self):
        cluster_id = None
        symbol = None
        annotation = None


def monomer_graph_analysis_check(input_file_name, monomer_file_name, output_prefix 
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
    print("Max head-to-tail distance: ", head_to_tail)
    print("Shortest read length: ", min_fasta_len)
    print("Clustering thresolds: ", thres_list)


    # IMPORT FASTA FILES #
    with open(input_fasta_file, 'r') as hin:
        for name, seq, qual in utils.readfq(hin):
            if len(seq) < min_fasta_len:
                print(">%s\n%s" % (name, seq), file = too_short_read)
                continue
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
            rng = rng.split("_")
            rng = int(rng[0]), int(rng[1])
            monomer_db.setdefault(rid, [])
            monomer_db[rid].append((rng, seq, orientation))

    print(len(seq_db), " sequences read." , "Reads with monomers:", len(monomer_db.keys()))


    # RUN OVER ALL READS THAT CONTAIN MONOMERS #
    for rid, monomers in monomer_db.items():
        aln_data = []
        range_list = []
        total_monomer_len = 0
        # Align all monomers on the read to each other.
        for i in range(len(monomers)):
            for j in range(i + 1, len(monomers)):
                t_seq = monomers[i][1]
                q_seq = monomers[j][1]
                aln = edlib.align(q_seq, t_seq)
                match_ratio = 1 - float(aln["editDistance"]) / min(len(q_seq), len(t_seq))  
                # Skip if no alignment
                if match_ratio >= 0.5:
                    # Store the alignment score for the monomer pair: (i, j, score)
                    aln_data.append((i, j, match_ratio))

            mono_range = monomers[i][0]
            # Add up app monomer ranges to calculate total monomer content in read
            total_monomer_len += (mono_range[1] - mono_range[0])
            # Store start-end coordinates of the monomer.
            range_list.append(mono_range[0])
            range_list.append(mono_range[1])

        is_regular = False
        inversion_detected = False
        missing_monomer = False


        # SCAN THE READ FOR AN HOR AT MULTIPLE CLUSTERING THRESHOLDS #
        for threshold in thres_list:
            G = nx.Graph()
            idt_in_clusters = []
            idt_out_clusters = []
            # Run over each monomer pair in aln_data to cluster monomers.
            for i, j, idt in aln_data:
                # Connect monomers i and j, if score is larger than threshold.
                # Monomers in the same cluster are treated as the same kind.
                if idt >= threshold:
                    G.add_edge(i, j)
                    idt_in_clusters.append(idt)
                else:
                    idt_out_clusters.append(idt)

            cluster_count = 0
            data = []
            data_c = {}
            l_seq = len(seq_db[rid])
            forward = False
            reverse = False

            # Check if monomers in either orientation exist in the read.
            for mono in monomers:
                orientation = mono[2]
                if orientation == "F":
                    forward = True
                else:
                    reverse = True
            # Mark read as INVERTED if monomers in both orientations exist.
            if forward and reverse:
                inversion_detected = True

            # Run over all clusters
            for C in nx.connected_components(G):
                # Run over all nodes (monomers) in cluster.
                for idx in C:
                    s, e = monomers[idx][0]
                    # Store (monomer_start, end, cluster_index) in 'data.'
                    # Store (monomer_start, cluster_index) in 'data_c.'
                    data.append((s, e, cluster_count))
                    data_c.setdefault(cluster_count, [])
                    data_c[cluster_count].append((s, cluster_count))
                cluster_count += 1

            if cluster_count <= 1 and not inversion_detected:
                # No or only one cluster. No HOR detected.\
                # Test with a different threshold.
                continue

            # Analyze clusters, if more than one detected.
            if cluster_count > 1:
                # Sort detected monomers by their coordinates.
                data.sort()
                # Create a symbolic HOR pattern, i.e., ABCDABCDABCD
                symbolic_pattern = ""
                for s, e, c in data:
                    symbolic_pattern += chr(65 + c)

                # Calculate head-to-tail distances between clustered monomers.
                x, e, y = zip(*data)
                x = np.array(x)
                head_to_head_intervals = x[1:] - x[:-1]
                head_to_tail_intervals = x[1:] - e[:-1] - 1

                # Calculate intervals between monomers in the same cluster
                c_intervals = []
                all_monomer_periods = []
                for c_index in data_c:
                    x = np.array([c[0] for c in data_c[c_index]])
                    x.sort()
                    # Periods: intervals between monomers of the same type
                    monomer_periods = x[1:] - x[:-1]
                    all_monomer_periods.extend(monomer_periods)

                min_monomer_period = min(all_monomer_periods)
                max_monomer_period = max(all_monomer_periods)
                median_monomer_period = round(np.median(all_monomer_periods))
                min_head_to_tail = min(head_to_tail_intervals)
                max_head_to_tail = max(head_to_tail_intervals)
                median_head_to_tail = round(np.median(head_to_tail_intervals))
                max_abs_head_to_tail = max(abs(head_to_tail_intervals))
                isolate_count = len(monomers) - len(data)
                HOR_start = min(range_list)
                HOR_end = max(range_list)
                monomeric_fraction_in_HOR = \
                    1.0*total_monomer_len/(HOR_end - HOR_start)

                # Normalize monomer periods by the median.
                if median_monomer_period > 0:
                    normalized_min_monomer_period = \
                        min_monomer_period/median_monomer_period
                    normalized_max_monomer_period = \
                        max_monomer_period/median_monomer_period
                else:
                    # Assign two out-of-range numbers if median is zero.
                    normalized_min_monomer_period = 0
                    normalized_max_monomer_period = 2

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

