#! /usr/bin/env python
import os
import subprocess
from . import utils

def parseHMMout(seq_db, input_file, monomer_file, orientation, min_monomer_len):

    ID = "NO ID"

    with open(input_file, 'r') as hin, open(monomer_file, 'w') as hout:
        for line in hin:
            l = line.rstrip('\n').split()    
            if l == []:
                continue
            if ">>" in line:      
                ID = l[1]
            else:
                if "!" in line:
                    low = int(l[12]) - 1
                    high = int(l[13]) - 1
                    if high - low + 1 >=  min_monomer_len:
                        if low != 0 and high < len(seq_db[ID]) - 1:
                            print(">%s/%d_%d/%s\n%s" % (ID, low + 1, high + 1, orientation, seq_db[ID][low:high+1]), file = hout)


def split_to_monomer_check(args):

    # args.min_monomer_len
    # monomers_file=in_seq_file.replace(".fa","_inferred_monomers.fa")

    # Call hmmsearch, build hmms based on consensus alignments
    if os.path.exists(args.output_file + ".hmmoutF.tbl"): 
        os.remove(args.output_file + ".hmmoutF.tbl")
    if os.path.exists(args.output_file + ".hmmoutF.out"): 
        os.remove(args.output_file + ".hmmoutF.out")
    subprocess.check_call(["hmmsearch", "--cpu", "8", 
        "--tblout", args.output_file + ".hmmoutF.tbl", 
        "-o", args.output_file + ".hmmoutF.out", 
        "--notextw", args.hmm_model_fwd, args.input_fasta_file])

    if os.path.exists(args.output_file + ".hmmoutR.tbl"): 
        os.remove(args.output_file + ".hmmoutR.tbl")
    if os.path.exists(args.output_file + ".hmmoutR.out"): 
        os.remove(args.output_file + ".hmmoutR.out")
    subprocess.check_call(["hmmsearch", "--cpu", "8", 
        "--tblout", args.output_file + ".hmmoutR.tbl",
        "-o", args.output_file + ".hmmoutR.out",
        "--notextw", args.hmm_model_rev, args.input_fasta_file])


    # import pdb; pdb.set_trace()
    seq_db = {}
    with open(args.input_fasta_file, 'r') as hin:
        for name, seq, qual in utils.readfq(hin):
            seq_db[name] = seq
    

    parseHMMout(seq_db, args.output_file + ".hmmoutF.out", 
        args.output_file + ".F.tmp", 
        'F', args.min_monomer_len)
    parseHMMout(seq_db, args.output_file + ".hmmoutR.out", 
        args.output_file + ".R.tmp", 
        'R', args.min_monomer_len)

    with open(args.output_file, 'w') as hout:
        subprocess.check_call(["cat", args.output_file + ".F.tmp", args.output_file + ".R.tmp"], stdout = hout)

    os.remove(args.output_file + ".F.tmp")
    os.remove(args.output_file + ".R.tmp")

    os.remove(args.output_file + ".hmmoutF.tbl")
    os.remove(args.output_file + ".hmmoutF.out")
    os.remove(args.output_file + ".hmmoutR.tbl")
    os.remove(args.output_file + ".hmmoutR.out")

