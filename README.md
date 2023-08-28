# DNA-LDPC-codes

This repository is the implementation of *Reducing Cost in DNA-based Data Storage by Sequence Analysis-aided Soft Information Decoding of Variable-Length Reads* (S.-J Park *et al.* (2023))


Thanks to libraries we used in this work..


* ***LDPC***: Library for encoding and decoding are used. It is available in https://github.com/radfordneal/LDPC-codes.  
* ***FLASH***: Library for merging paried-end reads from next-generation sequencing experiments. It is available in http://ccb.jhu.edu/software/FLASH/.  
* ***MUSCLE***: Library for multiple sequence alignment (MSA) algorithm. It is available in https://github.com/rcedgar/muscle.  

# Simulation environment
We run the DNA LDPC decoder with Anaconda on Windows 10 (https://www.anaconda.com/products/distribution).
Since .exe files are included, you have to run in Window version of anaconda.
```
            version
conda       4.12.0
python      3.8.8
```


# Codes

To run the decoder, go to the folder "ex_decoder" and use the following command.
This command runs 10 trials of decoding for random sampling number 72000 (From 72000_RS_0.txt to 72000_RS_9.txt).
```python
python decoder.py --rs 72000 --start 0 --end 10 --epsil 0.02
```

```
#Parameters

--rs:       Random sampling number
--start:    Iteration starting number of randomly sampled data
--end:      Iteration ending number of randomly sampled data
--epsil:    Initial epsilon value
```


# Required files for decoding

These files must be in the same folder to run the decoder. All files are included in the folder "ex_decoder".
```
MUSCLE.exe:                          file for multiple sequence alignment
decode_n18432_m2048_final.pchk:      pchk file for RS LDPC codes
decoder.py:                          main decoder for DNA-based data storage
def_func.py:                         several functions used for the decoding
ldpc.exe:                            LDPC decoder
pre_processing.py:                   file for index-based decoding
rs_dec.exe:                          file for Reed-Solomon code decoder
rs_dec_init.m:                       MATLAB file for 'rs_dec.exe'
72000_RS_0.txt:                      randomly sampled sequencing data
...
72000_RS_9.txt
72000_Q_0.txt:                       quality score for randomly sampled sequencing data
...
72000_Q_9.txt
codeword_n18432_m1860_1.txt:         original codewords
...
codeword_n18432_m1860_272.txt
```
The source code for creating "ldpc.exe" is in the folder "LDPC_dec".
You can obtain the file by running the folder "LDPC_dec" by the Release configuration setting, which creates an executable file and it is intended for infal, production-ready versions of your software.
Also, the source code for "rs_dec.exe" is in "rs_dec_init.m", and you can obtain the file by running "rs_dec_init.m" in the MATLAB.
