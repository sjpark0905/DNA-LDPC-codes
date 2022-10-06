"""
Copyright (c) 2022 by Seong-Joon Park
Coding and Cryptography Lab (CCL)
Department of Electrical and Computer Engineering,
Seoul National University, South Korea.
Email address: sjpark@ccl.snu.ac.kr
All Rights Reserved.
"""
# =============================================================================
# Perform LDPC decoding and re-decoding using LLR values
# =============================================================================

import os
import math
import numpy as np
from def_func import file_read, write_bat, write_codeword, re_index
from para import RS, re_th
# RS = int(input('Random Sampling Number: '))
# re_th = int(input('Re-decoding threshold(number): '))
# print('\n')


# def file_read(file,readtype):
#     with open(file+'.txt','r') as f:
#         data = f.readline()
#         data = data.split()
#         data = list(map(readtype,data))
    
#     return(data)
    
# def write_bat(file,soft,i):
#     with open(file, 'w') as f:
#     #bat = "ldpc 0.009000 0 0 0 4250 200 1 codeword%d LLR_inj%d bit_injection 0 1 0 0 1 768" % (i,i)
#         bat = "ldpc 0 0 0 7 200 1 codeword_n18432_m1860_%d %s decode_n18432_m2048_final 0 0 0 0" % (i,soft)
        
#         f.write(bat)
    
    
# def write_codeword(file,data):
#     with open(file+'.txt', 'w') as f:
#         for val in data:
#             f.write(val+' ')
    
    
# def re_index(x,th):
#     index = []
#     for i in range(18432):
#         if x[i]>th:
#             index.append(i+1)
#     return index
    
    
re_decode = [0 for i in range(18432)]
fail_DNA = []
# re_decode_index = []
# dict_error = {}

print('Start decoding')
for i in range(1,273):
    codeword_file = "codeword_n18432_m1860_%d" % i
    dec_file = "dec_codeword_n18432_m1860_%d" % i
    soft_input = "soft%d_n18432_m1860_%d" % (RS,i)
    write_bat("test.bat",soft_input,i)
    os.system('test.bat')
    
    dec_input = file_read(soft_input,float)
    dec_output = file_read(dec_file,int)
    codeword = file_read(codeword_file,int)
    v = 0
    
    for j in range(18432):
        
        if dec_input[j] >= 0:
            temp = 0
        else:
            temp = 1
        # print(temp)
        
        if dec_output[j] != temp:
            re_decode[j] += 1
            
        if dec_output[j] != codeword[j]:
            v += 1
        
            
    if v !=0:
        # print("틀렸어!!!") 
        # print(v)
        # dict_error[i] = v
        fail_DNA.append(i)
        print('index%d:   \terrors:%d' %(i,v))

print('First decoding finished\n')
# keyList = dict_error.keys()

# for item in keyList:
#     print('index: %d\terrors: %d' %(item,dict_error[item]))
re_decode_index = re_index(re_decode,re_th)  
print('Number of re-decoding strand: %d\n' % len(re_decode_index))
# print(re_decode_index)


for i in fail_DNA:
    max_llr = 0
    min_llr = 0
    soft_input = "soft%d_n18432_m1860_%d" % (RS,i)
    re_soft_input = "re_soft%d_n18432_m1860_%d" % (RS,i)
    dec_input = file_read(soft_input,float)
    
    for j in range(18432):
        if j in re_decode_index:
            dec_input[j] = 0.0
        elif dec_input[j]>=2*math.log((1-0.04)/0.04):
            max_llr+=1
            dec_input[j] = 30
        elif dec_input[j]<=-2*math.log((1-0.04)/0.04):
            min_llr+=1
            dec_input[j] = -30  
    print("# of max LLR:%d" % max_llr)
    print("# of min_LLR:%d\n" % min_llr)
    # for j in re_decode_index:
    #     dec_input[j] = 0.0
    
    dec_input = list(map(str,dec_input))
    write_codeword(re_soft_input, dec_input)
    

# =============================================================================
# Re-decode
# =============================================================================

# print(fail_DNA)

if fail_DNA:
    print("Start re-decoding")
fail_DNA2 = []
re_decode_num = 1

for i in fail_DNA:
    print("Re-decoding:     %d/%d" %(re_decode_num,len(fail_DNA)))
    re_decode_num += 1
    codeword_file = "codeword_n18432_m1860_%d" % i
    dec_file = "dec_codeword_n18432_m1860_%d" % i
    soft_input = "re_soft%d_n18432_m1860_%d" % (RS,i)
    write_bat("test.bat",soft_input,i)
    os.system('test.bat')
    
    dec_input = file_read(soft_input,float)
    dec_output = file_read(dec_file,int)
    codeword = file_read(codeword_file,int)
    v = 0
    
    for j in range(18432):
        
        if dec_input[j] >= 0:
            temp = 0
        else:
            temp = 1
        # print(temp)
        
        if dec_output[j] != temp:
            re_decode[j] += 1
            
        if dec_output[j] != codeword[j]:
            v += 1
        
            
    if v !=0:
        # print("틀렸어!!!") 
        # print(v)
        # dict_error[i] = v
        fail_DNA2.append(i)
        print('index:%d\terrors:%d' %(i,v))

print('Second decoding finished\n')
print('# ==============================================================================#')
print('#                    Results                                                    #')
print('# ==============================================================================#\n')
print("Random Sampling Number: %d" % RS)
print("Re-decoding threshold (number): %d\n" % re_th)
if not fail_DNA2:
    print("Decoding success")
    print("First decoding result:   %d/272" %(272-len(fail_DNA)))
    print("Second decoding result:  %d/272\n" %(272-len(fail_DNA2)))
    print("First decoding failure index")
    if not fail_DNA:
        print("None")
    else:
        print(fail_DNA)
    print("Second decoding failure index")
    if not fail_DNA2:
        print("None")
    else:
        print(fail_DNA2)
else:
    print("Decoding failure")
    print("First decoding result:   %d/272" %(272-len(fail_DNA)))
    print("Second decoding result:  %d/272\n" %(272-len(fail_DNA2)))
    print("First decoding failure index")
    if not fail_DNA:
        print("None")
    else:
        print(fail_DNA)
    print("Second decoding failure index")
    if not fail_DNA2:
        print("None")
    else:
        print(fail_DNA2)


        
    
    

# print(re_decode)

#dec_output[1] = 1
#print(dec_output[1])
#print(codeword[1])
#for i in range(len(codeword)):
#    if dec_output[i] != codeword[i]:
#        print("다르다!!")
#        break
    
    
    
    
    

