clear;
clc;

%%
A = sprintf('index.txt');
[FP] = fopen(A,'r');
ind = fscanf(FP,'%d');
ind = ind';

index_binary = zeros(length(ind)/32,32);

for i= 1:length(ind)/32
    
    index_binary(i,:) = ind(32*(i-1)+1:32*i);
    
end

gf_size = 4;

dec_index = zeros(size(index_binary,1),8);
for i=1:size(index_binary,1)
    a = zeros(1,8);
    for j=1:8

        temp = index_binary(i,gf_size*(j-1)+1:gf_size*j);
        a(j)=bin2dec(char(temp+48));
    end
    dec_index(i,:) = a;
end

dec_index = gf(dec_index,4);
[rs_dec,cnumerr] = rsdec(dec_index,8,4);

rs_dec = rs_dec.x;
dec_binary_index = zeros(size(index_binary,1),16);


for i = 1:size(index_binary,1)
    binary_temp = fliplr(de2bi(rs_dec(i,:)));
    binary_encoded = [];
    gf2bin = zeros(4,4);
    gf2bin(:,(4+1)-size(binary_temp,2):4)= binary_temp;


    for j=1:4
        binary_encoded = [binary_encoded gf2bin(j,:)];
    end

    dec_binary_index(i,:) = binary_encoded;
end

save dec_binary_index.mat dec_binary_index
save cnumerr.mat cnumerr
