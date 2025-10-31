function relation_matrix = hamming_matrix(data,delta)
    [~, attrinu] = size(data);
    hamming_distances = pdist(data, 'hamming')*attrinu; 
    hamming_matrix = squareform(hamming_distances);
    relation_matrix = exp(-hamming_matrix./delta); 
end