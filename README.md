# centr_utils
Utilities for centromere analysis based on [alpha-CENTAURI](https://github.com/volkansevim/alpha-CENTAURI).
The development of this repository was started to enhance my understanding on alpha satellite sequences.
But, I would like to extend it further.

## Dependency

### Binary 

HMMer

### Python package

networkx, edlib


## Install

```
https://github.com/friend1ws/centr_utils.git
cd centr_utils
python -m pip install . --user
```

## Command

### split_to_monomer
This is based on `chop_to_monomers.py` in alpha-CENTAURI.
```
centr_utils split_to_monomer [-h] [--min_monomer_len MIN_MONOMER_LEN]
  input_fasta_file output_file hmm_model_fwd hmm_model_rev
```

### monomer_graph_analysis
This is based on `monomer_graph_analysis.py` in alpha-CENTAURI.
```
centr_utils monomer_graph_analysis [-h]
  [--mean_monomer_len MEAN_MONOMER_LEN]
  [--head_to_tail_dist HEAD_TO_TAIL_DIST]
  [--min_fasta_len MIN_FASTA_LEN]
  [--thres_list THRES_LIST]
  input_fasta_file monomer_fasta_file output_prefix
```

