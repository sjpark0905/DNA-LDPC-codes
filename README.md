# DNA-LDPC-codes

This repository is the implementation of *Reducing Cost in DNA-based Data Storage by Sequence Analysis-aided Soft Information Decoding of Variable-Length Reads* (S.-J Park *et al.* (2023))


Thanks to libraries we used in this work..


* ***LDPC***: Library for encoding and decoding are used. It is available in https://github.com/radfordneal/LDPC-codes.  
* ***FLASH***: Library for merging paried-end reads from next-generation sequencing experiments. It is available in http://ccb.jhu.edu/software/FLASH/.  
* ***MUSCLE***: Library for multiple sequence alignment (MSA) algorithm. It is available in https://github.com/rcedgar/muscle.  

# Codes

```python
python decoder.py --rs 72000 --start 0 --end 10 --epsil 0.02
```
