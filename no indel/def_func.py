"""
Copyright (c) 2022 by Seong-Joon Park
Coding and Cryptography Lab (CCL)
Department of Electrical and Computer Engineering,
Seoul National University, South Korea.
Email address: sjpark@ccl.snu.ac.kr
All Rights Reserved.
"""

def file_read(file,readtype):
    all_data = []
    with open(file+'.txt','r') as f:
        while True:
            data = f.readline()
            
            if not data:
                break
            
            if readtype==str:
                all_data.append(data.strip('\n'))
            else:
                data = data.split()
                data = list(map(readtype,data))
                return(data)
            
    return(all_data)
    
def write_bat(file,soft,i):
    with open(file, 'w') as f:
    #bat = "ldpc 0.009000 0 0 0 4250 200 1 codeword%d LLR_inj%d bit_injection 0 1 0 0 1 768" % (i,i)
        bat = "ldpc 0 0 0 7 200 1 codeword_n18432_m1860_%d %s decode_n18432_m2048_final 0 0 0 0" % (i,soft)
        
        f.write(bat)
    
    
def write_codeword(file,data):
    with open(file+'.txt', 'w') as f:
        for val in data:
            f.write(val+' ')
    
    
def re_index(x,th):
    index = []
    for i in range(18432):
        # if x[i]>th and x[i]<160:
        if x[i]>th:
            #index.append(i+1)
            index.append(i)
    return index


class Fastq:
    def __init__(self, file):
        
        f = open(file+'.fastq','r')
        self.header = []
        self.seq = []
        self.qual = []
        
        for i,d in enumerate(f):
            if i%4 == 0:
                d = d.strip('\n')
                self.header.append(d)
            elif i%4 == 1:
                d = d.strip('\n')
                self.seq.append(d)
            elif i%4 == 3:
                d = d.strip('\n')
                self.qual.append(d)
            else:
                continue
            
def write_val(file,data):
    f = open(file+'.txt', 'w')
    for val in data:
        f.write(val+'\n')    
    f.close()
    
    
def DNA2binary(index):
    # b_index = []
    DNA_index_final = []
    for i in range(len(index)):
        b_index = ''
        for j in range(len(index[0])):
            if index[i][j] == 'A':
                temp = '0 0'
            elif index[i][j] == 'C':
                temp = '0 1'
            elif index[i][j] == 'G':
                temp = '1 0'
            else:
                temp = '1 1'
            b_index = b_index+temp+' '
            
        DNA_index_final.append(b_index)
        
    return DNA_index_final
        
        
def binary2decimal(x):
    temp = 0
    for i in range(len(x)):
        temp += x[i]*pow(2,len(x)-1-i)
    return temp