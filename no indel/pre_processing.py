"""
Copyright (c) 2022 by Seong-Joon Park
Coding and Cryptography Lab (CCL)
Department of Electrical and Computer Engineering,
Seoul National University, South Korea.
Email address: sjpark@ccl.snu.ac.kr or joonpark2247@gmail.com
All Rights Reserved.
"""

# =============================================================================
# Obtain the 18,432 true index values
# =============================================================================

import os
import math
import numpy as np
from def_func import binary2decimal

# def binary2decimal(x):
#     temp = 0
#     for i in range(len(x)):
#         temp += x[i]*pow(2,len(x)-1-i)
#     return temp
    
N = 18432
K = 16572
payload = 272

index_N = 16
index_K = 16
gf_size = 4
gf_index_N = index_N/gf_size
gf_index_K = index_K/gf_size

# index3 = np.zeros((index_K-2,pow(2,index_K-2)-1))
index = np.zeros((pow(2,index_K-2),index_K-2),int)                                  #0부터 2^14 까지를 binary로


# print(len(format(10,'b')))
# print(a)
for i in range(pow(2,index_K-2)):
    b_index = format(i,'b')
    index[i][index_K-2-len(b_index):index_K-2] = list(map(int,b_index))


DNA_index = np.zeros((pow(2,index_K-2),int((index_K-2)/2)),int)

for i in range(pow(2,index_K-2)):
    for j in range(int((index_K-2)/2)):
        temp = index[i][2*j:2*(j+1)]
        temp2 = temp[0]*2+temp[1]*1
        DNA_index[i][j] = temp2
                                                 
erase = []
DNA_index_final = np.empty((0,7),int)

for i in range(pow(2,index_K-2)):
    if DNA_index[i][2] == DNA_index[i][3] or DNA_index[i][5] == DNA_index[i][6]:        #3,4 번째, 6,7 번째 DNA가 다른 경우 HM을 만족하기 때문에.
        erase.append(i)
    else:
        DNA_index_final = np.append(DNA_index_final, [DNA_index[i]], axis=0)
        
binary_index = np.zeros((DNA_index_final.shape[0],index_K-2),int)

for i in range(DNA_index_final.shape[0]):
    for j in range(int((index_K-2)/2)):
        # b_index = format(DNA_index_final[i,j],'b')
        if DNA_index_final[i][j] == 0:
            temp = [0,0]
        elif DNA_index_final[i][j] == 1:
            temp = [0,1]
        elif DNA_index_final[i][j] == 2:
            temp = [1,0]
        else:
            temp = [1,1]
        binary_index[i][2*j:2*(j+1)] = temp



binary_index_final = np.zeros((2*DNA_index_final.shape[0],index_K),int)             #최종 18432개의 correct index

for i in range(DNA_index_final.shape[0]):
    for j in range(2):
        
        binary_index_final[2*i+j,0:index_K-2] = binary_index[i,:]
        binary_index_final[2*i+j][index_K-2:index_K] = [j,(sum(index[i])+j)%2]
        # binary_index_final[2*i+j,index_K-2] = (sum(index[i,:])+j)%2

decimal_index = []

for i in range(len(binary_index_final)):
    temp = 0
    for j in range(len(binary_index_final[0])):
        temp = binary2decimal(binary_index_final[i])
    decimal_index.append(temp) 
# binary_f!binary_final















