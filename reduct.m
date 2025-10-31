function select_attribute = reduct(data, delta, sigma)

[~, attrinu] = size(data);  

red = [];  
sig = 0;  
all_attri = 1:attrinu; 

for j=attrinu:-1:1
    remain_attri = setdiff(all_attri, red); 
    sig_tem = []; 
    for i = 1:length(remain_attri)
        red_i=union(red,remain_attri(i));
        red_i_matrix = hamming_matrix(data(:, red_i),delta); 
        temp = calculate_fuzzy_entropy(red_i_matrix);
        sig_tem = [sig_tem, temp];  
    end
    
    [x1, n1] = max(sig_tem);
    sig = [sig, x1];
    if abs(sig(end) - sig(end - 1)) > sigma
        red = [red, remain_attri(n1)];
    else
        break;
    end
end
    select_attribute=red;
end