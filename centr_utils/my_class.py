#! /usr/bin/env python

import statistics
import networkx as nx
import edlib

class Monomers(object):
    
    def __init__(self, origina_seq_id, original_seq_len):
        self.original_seq_id = origina_seq_id
        self.original_seq_len = original_seq_len

        self.monomer_list = [] 
        self.total_monomer_len = None

        self.aln_info = []

        self.is_regular = False
        self.inversion_detected = False
        self.missing_monomer = False

        self.cluster_count = 0
        self.cluster_list = []
        self.mean_identity_within_clusters = None
        self.mean_identity_between_clusters = None
        self.match_score_thres = None
        self.clustered_monomer_count = None
        self.clustered_monomer_len = None
        self.isolated_monomer_count = None

        self.min_head_to_tail = None
        self.max_head_to_tail = None
        self.median_head_to_tail = None
        self.max_abs_head_to_tail = None

        self.min_monomer_period = None
        self.max_monomer_period = None
        self.median_monomer_period = None
        self.normalized_min_monomer_period = None
        self.normalized_max_monomer_period = None


    def add_monomer(self, monomer):
        self.monomer_list.append(monomer)


    def comp_monomers(self):

        self.aln_data = []
        for i in range(len(self.monomer_list)):
            for j in range(i + 1, len(self.monomer_list)):
                monomer_i = self.monomer_list[i]
                monomer_j = self.monomer_list[j]
                aln = edlib.align(monomer_i.seq, momoner_j.seq)
                match_ratio = 1 - float(aln["editDistance"]) / min(len(monomer_i.seq), len(monomer_j.seq))
                # Skip if no alignment
                if match_ratio >= 0.5:
                    # Store the alignment score for the monomer pair: (i, j, score)
                    self.aln_data.append((i, j, match_ratio))

    
    def cluster_monomers(self, match_score_thres):

        G = nx.Graph()
        idt_in_clusters = []
        idt_out_clusters = []
        cluster_id = 0
        self.cluster_list = []
        clusterd_monomer_count = 0
        clustered_monomer_len = 0
        # self.cluster_count = 0
        # Run over each monomer pair in aln_data to cluster monomers.
        for i, j, idt in self.aln_data:
            # Connect monomers i and j, if score is larger than threshold.
            # Monomers in the same cluster are treated as the same kind.
            if idt >= match_score_thres:
                G.add_edge(i, j)
                idt_in_clusters.append(idt)
            else:
                idt_out_clusters.append(idt)

        # Run over all clusters
        for C in nx.connected_components(G):
            # Run over all nodes (monomers) in cluster.
            tcluster = Cluster(cluster_id)
            tcluster.symbol = chr(65 + cluster_id)
            self.cluster_list.append(tcluster_id)
            for idx in C:
                mono = self.monomer_list[idx]
                mono.cluster_id = cluster_id
                clustered_monomer_count = clustered_monomer_count + 1
                clustered_monomer_len = clustered_monomer_len + mono.end - mono.start
                # s, e = mono.start, mono.end
                # Store (monomer_start, end, cluster_index) in 'data.'
                # Store (monomer_start, cluster_index) in 'data_c.'
                # data.append((s, e, cluster_count))
                # data_c.setdefault(cluster_count, [])
                # data_c[cluster_count].append((s, cluster_count))
                cluster_id = cluster_id + 1

        self.cluster_count = cluster_id
        self.mean_identity_within_clusters = statistics.mean(idt_in_clusters)
        self.mean_identity_between_clusters = statistics.mean(idt_out_clusters)

        self.clustered_monomer_count = clustered_monomer_count
        self.clustered_monomer_len = clustered_monomer_len
        self.isolated_monomer_count = len(self.monomer_list) - self.clustered_monomer_count
        self.match_score_thres = match_score_thres
        

    def inversion_check(self):

        is_forward = False
        is_reverse = False
        # Check if monomers in either orientation exist in the read.
        for monomer in self.monomer_list:
            if monomer.orientation == 'F':
                is_forward = True
            else:
                is_reverse = True
        # Mark read as INVERTED if monomers in both orientations exist.
        if is_forward and is_reverse:
            self.inversion_detected = True


    def set_head_to_tail_stat(self):

        # assume monomer_list are sorted by the coordinate of start position
        start_list_clustered = []
        end_list_clustered = []
        for mono in monomer_list:
            if mono.clustered_id is not None
                start_list_clustered.append(mono.start)
                end_list_clustered.append(mono.end)

        head_to_tail_intervals = start_list_clustered[1:] - end_list_clustered[:-1] - 1
   
        self.min_head_to_tail = min(head_to_tail_intervals)
        self.max_head_to_tail = max(head_to_tail_intervals)
        self.median_head_to_tail = statistics.median(head_to_tail_intervals)
        self.max_abs_head_to_tail = max(abs(head_to_tail_intervals))

    
    def set_monomer_period_stat(self):

        # assume monomer_list are sorted by the coordinate of start position
        cid2start = {}
        for mono in monomer_list:
            if mono.clustered_id is not None:
                if cid not in cid2start: cid2start = []
                cid2start[mono.clustered_id].append(mono.start)

        all_monomer_periods = []
        for cid in sorted(cid2start):
            monomer_periods = cid2start[cid][1:] - cid2start[cid][:-1]
            all_monomer_periods.extend(monomer_periods)

        self.min_monomer_period = min(all_monomer_periods)
        self.max_monomer_period = max(all_monomer_periods)
        self.median_monomer_period = statistics.median(all_monomer_periods)

        # Normalize monomer periods by the median.
        if self.median_monomer_period > 0:
            self.normalized_min_monomer_period = \
                self.min_monomer_period / self.median_monomer_period
            self.normalized_max_monomer_period = \
                self.max_monomer_period / self.median_monomer_period
        else:
            # Assign two out-of-range numbers if median is zero.
            self.normalized_min_monomer_period = 0
            self.normalized_max_monomer_period = 2

        """
        HOR_start = min(range_list)
        HOR_end = max(range_list)
        monomeric_fraction_in_HOR = \
            1.0*total_monomer_len/(HOR_end - HOR_start)
        """


class Monomer(object):
    
    def __init__(self, start, end, seq, orientation):
        self.start = start
        self.end = end
        self.seq = seq
        self.orientation = orientation
        self.cluster_id = None

class Cluster(ojbect):

    def __init__(self, cluster_id):
        self.cluster_id = cluster_id
        self.symbol = None
        self.annotation = None
