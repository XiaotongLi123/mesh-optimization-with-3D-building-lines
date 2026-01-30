function [neigh_idx]=getNeighKidx(idx,this_idx,neigh_k)
    idx_num=size(idx,1);
    idx=[idx;idx;idx];
    neigh_idx=idx(idx_num+this_idx-neigh_k:idx_num+this_idx+neigh_k,1);
end