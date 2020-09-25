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


def split_to_monomer_check(input_fasta_file, output_file, hmm_model_fwd, hmm_model_rev, min_monomer_len):


    # Call hmmsearch, build hmms based on consensus alignments
    if os.path.exists(output_file + ".hmmoutF.tbl"): 
        os.remove(output_file + ".hmmoutF.tbl")
    if os.path.exists(output_file + ".hmmoutF.out"): 
        os.remove(output_file + ".hmmoutF.out")
    subprocess.check_call(["hmmsearch", "--cpu", "8", 
        "--tblout", output_file + ".hmmoutF.tbl", 
        "-o", output_file + ".hmmoutF.out", 
        "--notextw", hmm_model_fwd, input_fasta_file])

    if os.path.exists(output_file + ".hmmoutR.tbl"): 
        os.remove(output_file + ".hmmoutR.tbl")
    if os.path.exists(output_file + ".hmmoutR.out"): 
        os.remove(output_file + ".hmmoutR.out")
    subprocess.check_call(["hmmsearch", "--cpu", "8", 
        "--tblout", output_file + ".hmmoutR.tbl",
        "-o", output_file + ".hmmoutR.out",
        "--notextw", hmm_model_rev, input_fasta_file])


    # import pdb; pdb.set_trace()
    seq_db = {}
    with open(input_fasta_file, 'r') as hin:
        for name, seq, qual in utils.readfq(hin):
            seq_db[name] = seq
    

    parseHMMout(seq_db, output_file + ".hmmoutF.out", 
        output_file + ".F.tmp", 
        'F', min_monomer_len)
    parseHMMout(seq_db, output_file + ".hmmoutR.out", 
        output_file + ".R.tmp", 
        'R', min_monomer_len)

    with open(output_file, 'w') as hout:
        subprocess.check_call(["cat", output_file + ".F.tmp", output_file + ".R.tmp"], stdout = hout)

    os.remove(output_file + ".F.tmp")
    os.remove(output_file + ".R.tmp")

    os.remove(output_file + ".hmmoutF.tbl")
    os.remove(output_file + ".hmmoutF.out")
    os.remove(output_file + ".hmmoutR.tbl")
    os.remove(output_file + ".hmmoutR.out")

