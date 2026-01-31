function [result] = getQuantile(data, q)
    data_sort = sort(data, 'ascend');
    idx = ceil(size(data, 1) * q);
    result = data_sort(idx, 1);
end