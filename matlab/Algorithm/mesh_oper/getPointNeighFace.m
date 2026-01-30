function [neigh_face_idx]=getPointNeighFace(point,PFneighbor)
    % Find all faces that are adjacent to a given set of points
    point_num=size(point,1);
    neigh_face_idx_ori=zeros(19*point_num,1);
    neigh_face_num=0;
    for i=1:point_num
        for j=1:PFneighbor(point(i,1),1)
            if ~any(ismember(neigh_face_idx_ori(1:neigh_face_num,1),PFneighbor(point(i,1),j+1)))
                neigh_face_num=neigh_face_num+1;
                neigh_face_idx_ori(neigh_face_num,1)=PFneighbor(point(i,1),j+1);
            end
        end
    end
    neigh_face_idx=neigh_face_idx_ori(1:neigh_face_num,1);
end