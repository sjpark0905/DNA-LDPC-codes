"""
Copyright (c) 2022 by Seong-Joon Park
Coding and Cryptography Lab (CCL)
Department of Electrical and Computer Engineering,
Seoul National University, South Korea.
Email address: joon2247@snu.ac.kr;joonpark2247@gmail.com
All Rights Reserved.
"""

import math
import numpy as np
import os
from scipy import io
from pre_processing import binary2decimal, decimal_index
from def_func import write_val, DNA2binary, Fastq, write_codeword, file_read, write_bat, re_index, edit_dist
import argparse
import time
import random


parser = argparse.ArgumentParser(description='Decoding of the sequenced DNA data')

parser.add_argument('--rs',          type=int,   default=70000, help='Random sampling number')
parser.add_argument('--start',     type=int,   default=0, help='Iteration start number')
parser.add_argument('--end',     type=int,   default=10, help='Iteration end number')
parser.add_argument('--epsil',     type=float,   default=0.03, help='Epsilon value')

args = parser.parse_args()


print(args)

RS = args.rs
start = args.start
end = args.end
epsil = args.epsil


with open('alignment.bat', 'w') as f:
    bat = "MUSCLE -align align_input.txt -output align_output.txt"
    f.write(bat)
    

for iter_num2 in range(start,end):
    index = list()
    RS_seq = list()
    time_RS = time.time()
    RS_file = "%d_RS_%d" %(RS,iter_num2)
    RS_file_qual = "%d_RS_Q_%d" %(RS,iter_num2)
    
    if os.path.exists(RS_file+".txt"):
        print("************** Read random sampling file! **************")
        RS_seq = file_read(RS_file,str)
        RS_qual = file_read(RS_file_qual,str)
    else:
        print("************** No random sampling file! **************")
        break
    
    for i in range(len(RS_seq)):
        index.append(RS_seq[i][0:16])
    
    time_RS_dec = time.time()
    index_binary = DNA2binary(index)
    write_val('index',index_binary)
    
    
    with open('RS_decoder.bat', 'w') as f:
        bat = "rs_dec"
        f.write(bat)
    os.system('RS_decoder.bat')
    
    print('RS decoding finished')
    print("RS decoding time:", time.time()-time_RS_dec,"sec")
   
    
    mat_file1 = io.loadmat('dec_binary_index.mat')
    RS_dec_index = mat_file1['dec_binary_index']
    
    mat_file2 = io.loadmat('cnumerr.mat')
    cnumerr = mat_file2['cnumerr']
    
    RS_dec_index_final = np.empty((0,16),int)
    RS_dec_codeword_final = np.empty((0,len(RS_seq[0])),int)
    RS_qual_final = list()
    
    for i in range(len(RS_seq)):
        if cnumerr[i][0] == 0 or cnumerr[i][0] == 1 or cnumerr[i][0] == 2:
            RS_dec_index_final = np.append(RS_dec_index_final, [RS_dec_index[i]], axis=0)
            RS_dec_codeword_final = np.append(RS_dec_codeword_final, [RS_seq[i]])
            RS_qual_final.append(ord(RS_qual[i]))
        else:
            continue
    
    decimal_index_cand = []
    
    for i in range(len(RS_dec_index_final)):
        temp = 0
        for j in range(len(RS_dec_index_final[0])):
            temp = binary2decimal(RS_dec_index_final[i])
        decimal_index_cand.append(temp)
            
    without_error_index2 = []    
    index_codeword = []                                                 
    index_DNA = []                                                
    error_index =[]                                                    
                                                        
    
    

    for i,index_cand in enumerate(decimal_index_cand):
        if index_cand not in decimal_index:
            error_index.append(index_cand)                              
        else:
            without_error_index2.append(index_cand)
            index_DNA.append((decimal_index_cand[i],RS_dec_codeword_final[i][16:len(RS_dec_codeword_final[i])],RS_qual_final[i]))
            
    index_DNA = sorted(index_DNA, key=lambda x: x[0])
    index_codeword = sorted(index_codeword, key=lambda x: x[0])         
    without_error_index2 = sorted(without_error_index2) 
    time_LLR_cal = time.time()
    erase_index = []
    DNA_cand = []
    LLR_cand =[]
    q_272 = []
    DNA_LLR = [0 for i in range(272)]
    DNA_LR = [0 for i in range(272)]
    DNA_index_LLR_final = []
    DNA_index_LR_final = []
    aligned_DNA=list()
    count_0 = 0
    count_1 = 0

    q_count_0 = 0
    q_count_1 = 0

    v1 = 0          
    v2 = 0          
    num_test = 0
    index_test2 = []
    error_align = []
    error_align2 =[]
    error_temp = []

    without_error_index=list()

    print('************** Start MSA and LLR calculation **************\n')

    while v1 <= len(index_DNA):
        
        if decimal_index[v2] not in without_error_index2:
            erase_index.append(decimal_index[v2])
            v2 += 1
            continue

        if v1 == len(index_DNA):       
        
            r_q_272 = list()
            r_DNA_cand = list()
            q_272_2 = list()
        
            temp_order = list()
            aligned_DNA = list()
            if len(DNA_cand)!=1:
                
                test_len = 0
                for num_1 in range(len(DNA_cand)):
                    if len(DNA_cand[num_1]) == 136:
                        test_len +=1
                if test_len == len(DNA_cand):
                    num_test+=1
                    r_DNA_cand = DNA_cand
                    r_q_272 = q_272
                    LLR_cand = DNA2binary(r_DNA_cand)
                    
                else:
                        
                        
                    same_seq = list()
                    for i in range(len(DNA_cand)):
                        for k in range(i+1,len(DNA_cand)):
                            temp = edit_dist(DNA_cand[i], DNA_cand[k])
                            if temp < 15:
                                same_seq.append(i)
                                same_seq.append(k)
                    for i in range(len(np.unique(same_seq))):
                        r_DNA_cand.append(DNA_cand[np.unique(same_seq)[i]])
                        q_272_2.append(q_272[np.unique(same_seq)[i]])
                    if len(r_DNA_cand)==0:
                        align_DNA = []
                        aligned_DNA =[]
                        DNA_cand = []
                        LLR_cand =[]
                        q_272 = []
                        DNA_LLR = [0 for i in range(272)]
                        DNA_LR = [0 for i in range(272)]
                        v2 += 1
                        continue
                    
                    
                    
                    with open('align_input.txt', 'w') as f:
                        for val_num,val in enumerate(r_DNA_cand):
                            write_num = ">%d\n%s\n" %(val_num,val)
                            f.write(write_num)
    
                    os.system('alignment.bat')
                    align_DNA = file_read('align_output',str)
    
                    temp =''
                    error_q = list()
                    for num_val, val in enumerate(align_DNA):
                        
                        if num_val%3==0:
                            temp_order = val[1:]
                            
                            continue
                        elif num_val%3==1:
                            temp +=val
    
                        elif num_val%3==2:
                            temp +=val
    
                            if len(temp)!=136:
                                error_q.append([q_272_2[int(temp_order)],temp[len(temp)-1]])
                                temp =''
                                continue
                            r_q_272.append(q_272_2[int(temp_order)])
                            aligned_DNA.append(temp)
                            temp_order = list()
                            temp =''
                        if len(DNA_cand) != len(aligned_DNA):
                            error_temp.append(aligned_DNA)
                    LLR_cand = DNA2binary(aligned_DNA)
            else:
                r_DNA_cand = DNA_cand
                r_q_272 = q_272
                if len(r_DNA_cand[0])<136:
                    LLR_cand = DNA2binary(r_DNA_cand)
                    LLR_cand[0]=LLR_cand[0].replace(" ","")
                    
                    if r_q_272[0] >63:
                        if LLR_cand[0][len(LLR_cand[0])-1] == '0':
                            DNA_LLR[271] = math.log((1-epsil)/epsil)
                        else:
                            DNA_LLR[271] = -math.log((1-epsil)/epsil)
                    
                    DNA_LR = np.exp(DNA_LLR)
                    without_error_index.append(decimal_index[v2])
                    DNA_index_LLR_final.append((decimal_index[v2],DNA_LLR))
                    DNA_index_LR_final.append((decimal_index[v2],DNA_LR))
                    
                    align_DNA = []
                    aligned_DNA =[]
                    DNA_cand = []
                    LLR_cand =[]
                    q_272 = []
                    DNA_LLR = [0 for i in range(272)]
                    DNA_LR = [0 for i in range(272)]
                    
                    v2 += 1
                    continue
                else:
                    
                    LLR_cand = DNA2binary(r_DNA_cand)

            for i in range(272):
                
                if len(LLR_cand)==0:
                    for error_len in range(len(error_q)):
                        temp2 =''
                        if error_q[error_len][0]>63:
                            temp2 = DNA2binary(error_q[error_len][1])
                            temp2 = temp2[0].replace(" ","")
                            
                            if temp2[1] == '0':
                                count_0 +=1
                            else:
                                count_1 +=1
                    DNA_LLR[271] = (count_0-count_1)*math.log((1-epsil)/epsil)
                    
                    
                    
                    DNA_LR[i] = np.exp(DNA_LLR[i])
                    q_count_0 = 0
                    q_count_1 = 0
                    count_0 = 0
                    count_1 = 0
                    error_q = list()
                    break
                
                
                for j in range(len(LLR_cand)):
                    LLR_cand[j]=LLR_cand[j].replace(" ","")
                    if (i==271) and (r_q_272[j] < 53):
                        continue
                
                   
                    if LLR_cand[j][i] == '0':
                        count_0 += 1
                        q_count_0 += r_q_272[j]
                    else:
                        count_1 += 1
                        q_count_1 += r_q_272[j]

                if (i==271) and (count_0==1) and (count_1==1):
                    if (q_count_0 < 53) and (q_count_1 >=63):
                        DNA_LLR[i] = -2*math.log((1-epsil)/epsil)
                    elif (q_count_0 >=63) and (q_count_1 < 53):
                        DNA_LLR[i] = 2*math.log((1-epsil)/epsil)
                    else:
                        DNA_LLR[i] = 0

                else:
                    DNA_LLR[i] = (count_0-count_1)*math.log((1-epsil)/epsil) 

                DNA_LR[i] = np.exp(DNA_LLR[i])
                q_count_0 = 0
                q_count_1 = 0
                count_0 = 0
                count_1 = 0
            without_error_index.append(decimal_index[v2])
            DNA_index_LLR_final.append((decimal_index[v2],DNA_LLR))
            DNA_index_LR_final.append((decimal_index[v2],DNA_LR))
            break

        if (decimal_index[v2]==index_DNA[v1][0]):
            DNA_cand.append(index_DNA[v1][1])
            q_272.append(int(index_DNA[v1][2]))
            v1 += 1
        else:
            
            r_q_272 = list()
            r_DNA_cand = list()
            q_272_2 = list()
        
            temp_order = list()
            aligned_DNA = list()
            if len(DNA_cand)!=1:
                
                test_len = 0
                for num_1 in range(len(DNA_cand)):
                    if len(DNA_cand[num_1]) == 136:
                        test_len +=1
                if test_len == len(DNA_cand):
                    num_test+=1
                    r_DNA_cand = DNA_cand
                    r_q_272 = q_272
                    LLR_cand = DNA2binary(r_DNA_cand)
                    
                else:
                        
                        
                    same_seq = list()
                    for i in range(len(DNA_cand)):
                        for k in range(i+1,len(DNA_cand)):
                            temp = edit_dist(DNA_cand[i], DNA_cand[k])
                            if temp < 15:
                                same_seq.append(i)
                                same_seq.append(k)
                    for i in range(len(np.unique(same_seq))):
                        r_DNA_cand.append(DNA_cand[np.unique(same_seq)[i]])
                        q_272_2.append(q_272[np.unique(same_seq)[i]])
                    if len(r_DNA_cand)==0:
                        align_DNA = []
                        aligned_DNA =[]
                        DNA_cand = []
                        LLR_cand =[]
                        q_272 = []
                        DNA_LLR = [0 for i in range(272)]
                        DNA_LR = [0 for i in range(272)]
                        v2 += 1
                        continue
                    
                    
                    
                    with open('align_input.txt', 'w') as f:
                        for val_num,val in enumerate(r_DNA_cand):
                            write_num = ">%d\n%s\n" %(val_num,val)
                            f.write(write_num)
    
                    os.system('alignment.bat')
                    align_DNA = file_read('align_output',str)
    
                    temp =''
                    error_q = list()
                    for num_val, val in enumerate(align_DNA):
                        if num_val%3==0:
                            temp_order = val[1:]
                            
                            continue
                        elif num_val%3==1:
                            temp +=val
    
                        elif num_val%3==2:
                            temp +=val
    
                            if len(temp)!=136:
                                error_q.append([q_272_2[int(temp_order)],temp[len(temp)-1]])
                                temp =''
                                
                                continue
                            r_q_272.append(q_272_2[int(temp_order)])
                            aligned_DNA.append(temp)
                            temp_order = list()
                            temp =''
                        if len(DNA_cand) != len(aligned_DNA):
                            error_temp.append(aligned_DNA)
                    LLR_cand = DNA2binary(aligned_DNA)
                
            else:
                r_DNA_cand = DNA_cand
                r_q_272 = q_272
                if len(r_DNA_cand[0])<136:
                    LLR_cand = DNA2binary(r_DNA_cand)
                    LLR_cand[0]=LLR_cand[0].replace(" ","")
                    
                    if r_q_272[0] >63:
                        if LLR_cand[0][len(LLR_cand[0])-1] == '0':
                            DNA_LLR[271] = math.log((1-epsil)/epsil)
                        else:
                            DNA_LLR[271] = -math.log((1-epsil)/epsil)
                    
                    DNA_LR = np.exp(DNA_LLR)
                    without_error_index.append(decimal_index[v2])
                    DNA_index_LLR_final.append((decimal_index[v2],DNA_LLR))
                    DNA_index_LR_final.append((decimal_index[v2],DNA_LR))
                    
                    align_DNA = []
                    aligned_DNA =[]
                    DNA_cand = []
                    LLR_cand =[]
                    q_272 = []
                    DNA_LLR = [0 for i in range(272)]
                    DNA_LR = [0 for i in range(272)]
                    
                    v2 += 1
                    continue
                else:
                    
                    LLR_cand = DNA2binary(r_DNA_cand)

            for i in range(272):
                
                if len(LLR_cand)==0:
                    for error_len in range(len(error_q)):
                        temp2 =''
                        if error_q[error_len][0]>63:
                            temp2 = DNA2binary(error_q[error_len][1])
                            temp2 = temp2[0].replace(" ","")
                            
                            if temp2[1] == '0':
                                count_0 +=1
                            else:
                                count_1 +=1
                    DNA_LLR[271] = (count_0-count_1)*math.log((1-epsil)/epsil)
                    
                    
                    
                    DNA_LR[i] = np.exp(DNA_LLR[i])
                    q_count_0 = 0
                    q_count_1 = 0
                    count_0 = 0
                    count_1 = 0
                    error_q = list()
                    break
                
                
                for j in range(len(LLR_cand)):
                    LLR_cand[j]=LLR_cand[j].replace(" ","")
                    if (i==271) and (r_q_272[j] < 53):
                        continue
                
                   
                    if LLR_cand[j][i] == '0':
                        count_0 += 1
                        q_count_0 += r_q_272[j]
                    else:
                        count_1 += 1
                        q_count_1 += r_q_272[j]

                if (i==271) and (count_0==1) and (count_1==1):

                    if (q_count_0 < 53) and (q_count_1 >=63):
                        DNA_LLR[i] = -2*math.log((1-epsil)/epsil)
                    elif (q_count_0 >=63) and (q_count_1 < 53):
                        DNA_LLR[i] = 2*math.log((1-epsil)/epsil)
                    else:
                        DNA_LLR[i] = 0

                else:
                    DNA_LLR[i] = (count_0-count_1)*math.log((1-epsil)/epsil)  

                DNA_LR[i] = np.exp(DNA_LLR[i])
                q_count_0 = 0
                q_count_1 = 0
                count_0 = 0
                count_1 = 0
            
            without_error_index.append(decimal_index[v2])
            DNA_index_LLR_final.append((decimal_index[v2],DNA_LLR))
            DNA_index_LR_final.append((decimal_index[v2],DNA_LR))
            align_DNA = []
            aligned_DNA =[]
            r_DNA_cand = []
            DNA_cand = []
            LLR_cand =[]
            q_272 = []
            DNA_LLR = [0 for i in range(272)]
            DNA_LR = [0 for i in range(272)]
            v2 += 1
    print('************** LLR calculation ended **************\n')


    for index_cand in decimal_index:
        if index_cand not in np.unique(without_error_index):
            DNA_index_LLR_final.append((index_cand, [0 for i in range(272)]))   
            DNA_index_LR_final.append((index_cand, [1 for i in range(272)]))      
    DNA_index_LLR_final = sorted(DNA_index_LLR_final, key=lambda x: x[0])         
    DNA_index_LR_final = sorted(DNA_index_LR_final, key=lambda x: x[0])        

    DNA_LLR_final =[]
    DNA_LR_final =[]

    for i in range(18432):
        DNA_LLR_final.append(DNA_index_LLR_final[i][1])
        DNA_LR_final.append(DNA_index_LR_final[i][1])

    LDPC_LLR_input = [0 for i in range(18432)]
    
    for i in range(272):
        for j in range(18432):
            LDPC_LLR_input[j] = DNA_LLR_final[j][i]
        LDPC_LLR_input = list(map(str,LDPC_LLR_input))
        soft_output = "soft%d_n18432_m1860_%d" % (RS,i+1)
        write_codeword(soft_output, LDPC_LLR_input)
        
    print('LLR calculation ended\n')
    print("LLR calculation time:", time.time()-time_LLR_cal,"sec")
    
    
    time_first_dec = time.time()
    
    
    re_decode = [0 for i in range(18432)]
    fail_DNA = []
    re_decode_index = []
    
    # ===============================================================================================
    # First decoding start
    # ===============================================================================================
    
    print('************** Start LDPC decoding **************\n')
    for i in range(1,273):
        codeword_file = "codeword_n18432_m1860_%d" % i
        dec_file = "dec_codeword_n18432_m1860_%d" % i
        soft_input = "soft%d_n18432_m1860_%d" % (RS,i)
        write_bat("soft_decoder.bat",soft_input,i)
        os.system('soft_decoder.bat')
        
        dec_input = file_read(soft_input,float)
        dec_output = file_read(dec_file,int)
        codeword = file_read(codeword_file,int)
        v = 0
        
        for j in range(18432):
            
            if dec_input[j] >= 0:
                temp = 0
            else:
                temp = 1
            
            if dec_output[j] != temp:
                re_decode[j] += 1
                
            if dec_output[j] != codeword[j]:
                v += 1
            
                
        if v !=0:
            fail_DNA.append(i)
            print('index: %d   \terrors: %d' %(i,v))
    
    print('First decoding finished\n')
    print("First decoding time:", time.time()-time_first_dec,"sec")
    
    # ===============================================================================================
    # Second decoding calculation
    # ===============================================================================================
    error_iter = list()
    iter_sec_dec = 0
    erasure_index = re_index(re_decode,140)  
    fail_DNA2 = fail_DNA
    
    epsil2 = epsil-0.0005
    while fail_DNA2 and epsil2 > 0.001:
        num_fail_DNA = len(fail_DNA2)
        iter_sec_dec = iter_sec_dec+1
        print("epsilon", epsil2)
        print(math.log((1-epsil2)/epsil2))
        for i in fail_DNA2:
            soft_input = "soft%d_n18432_m1860_%d" % (RS,i)
            re_soft_input = "re_soft%d_n18432_m1860_%d" % (RS,i)
            dec_input = file_read(soft_input,float)
                    
            for j in range(18432):
                 if dec_input[j] == 0:
                     continue
                 else:
                     dec_input[j] = dec_input[j]* math.log((1-epsil2+0.0005)/(epsil2-0.0005))/math.log((1-epsil)/epsil)
                     
        
            dec_input = list(map(str,dec_input))
            write_codeword(re_soft_input, dec_input)
        epsil2 = epsil2-0.0005
        
        # ===============================================================================================
        # Second decoding start
        # ===============================================================================================
        

        re_decode = [0 for i in range(18432)]
        
        if fail_DNA2:
            print("************** Start LDPC re-decoding **************\n")
            print('Number of erasure strand: %d\n' % len(erasure_index))

        re_decode_num = 0
        time_sec_dec = time.time()
        
        for i in fail_DNA2:
            re_decode_num += 1
            codeword_file = "codeword_n18432_m1860_%d" % i
            dec_file = "dec_codeword_n18432_m1860_%d" % i
            soft_input = "re_soft%d_n18432_m1860_%d" % (RS,i)
            write_bat("soft_decoder.bat",soft_input,i)
            os.system('soft_decoder.bat')
            
            dec_input = file_read(soft_input,float)
            dec_output = file_read(dec_file,int)
            codeword = file_read(codeword_file,int)
            v = 0
            
                    
            for j in range(18432):
                
                if dec_input[j] >= 0:
                    temp = 0
                else:
                    temp = 1
                if dec_output[j] != codeword[j]:
                    v += 1
            if v !=0:
                for j in range(18432):
                    if dec_output[j] != temp:
                        re_decode[j] += 1
                    
                
                    
            print("%d/%d\t" %(re_decode_num,num_fail_DNA), "index: %d\terrors: %d\n" %(i,v),)
            fail_DNA2 = []
            if v !=0:
                fail_DNA2.append(i)
        print('Second decoding finished')
        print("Second decoding time:", time.time()-time_sec_dec,"sec\n")
    # ===============================================================================================
    # Third decoding calculation
    # ===============================================================================================
    num_fail_DNA = len(fail_DNA2)
    
    
    if not fail_DNA2:
        f = open('o_%d_%d_%f_result.txt' %(RS,iter_num2,epsil), 'w')
    else:
        f = open('x_%d_%d_%f_result.txt' %(RS,iter_num2,epsil), 'w')
    
    f.write('==============================================================================\n')
    f.write('                               Results                                        \n')
    f.write('==============================================================================\n')
    f_time = 'Total time: %f sec\n' %(time.time()-time_RS)
    f.write(f_time)
    f_RS = "Random Sampling Number: %d\n" % RS
    f.write(f_RS)
    if not fail_DNA2:
        f.write("Decoding success\n\n")
        f1 = "First decoding result:   %d/272\n" %(272-len(fail_DNA))
        f.write(f1)
        f2 = "Second decoding result:  %d/272\n" %(272-len(fail_DNA2))
        f.write(f2)
        f22 = "Second decoding iteration number:  %d\n" % iter_sec_dec
        f.write(f22)
        
        f.write("First decoding failure index: ")
        if not fail_DNA:
            f.write("None\n")
        else:
            for val in fail_DNA:
                f.write(str(val)+' ')
            f.write('\n')
        f.write("Second decoding failure index: ")
        if not fail_DNA2:
            f.write("None\n")
        else:
            for val in fail_DNA2:
                f.write(str(val)+' ')
            f.write('\n')
        
    else:
        f.write("Decoding failure\n\n")
        f4 = "First decoding result:\t%d/272\n" %(272-len(fail_DNA))
        f.write(f4)
        f5 = "Second decoding result:\t%d/272\n" %(272-len(fail_DNA2))
        f.write(f5)
        
        f.write("First decoding failure index: ")
        if not fail_DNA:
            f.write("None\n")
        else:
            for val in fail_DNA:
                f.write(str(val)+' ')
            f.write('\n')
        f.write("Second decoding failure index: ")
        if not fail_DNA2:
            f.write("None\n")
        else:
            for val in fail_DNA2:
                f.write(str(val)+' ')
            f.write('\n')
        
        
    f.close()
    
    
