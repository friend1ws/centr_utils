#! /usr/bin/env python3

import sys, re
from centr_utils.utils import readfq, reverse_complement

input_file = sys.argv[1]
hor_list = sys.argv[2]

gid2alpha = {}
with open(hor_list, 'r') as hin:
    for line in hin:
        F = line.rstrip('\n').split(' ')
        gid2alpha[F[1]] = F[0]


def flush_info(temp_gid, temp_seq, temp_fr):

    temp_fr_2 = list(set(temp_fr))
    if len(temp_fr_2) > 1:
        print("Inconsistent orientation!", file = sys.stderr)
        sys.exit(1)
    if temp_fr_2[0] == 'R': temp_seq = reversed(temp_seq)

    temp_gid = re.sub(r'.\d+$', '', temp_gid)
    alpha_id = gid2alpha[temp_gid]
    tind = 0
    for tseq in temp_seq:
        if temp_fr_2[0] == 'R': tseq = reverse_complement(tseq)
        print(">%s_%s_%d\n%s" % (alpha_id, temp_gid, tind, tseq))
        tind = tind + 1

temp_gid = ''
temp_seq = []
temp_fr = []

with open(input_file, 'r') as hin:
    for name, seq, qual in readfq(hin):
    
        gid, reg, fr = name.split('/')
        if gid != temp_gid:
            if temp_gid != '':

                flush_info(temp_gid, temp_seq, temp_fr)

            temp_gid = gid
            temp_seq = []
            temp_fr = []


        temp_seq.append(seq)
        temp_fr.append(fr)


flush_info(temp_gid, temp_seq, temp_fr)
