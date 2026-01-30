function [neigh_face_id]=findPointNeighFace(face,point)
    % Find the triangular faces adjacent to the current point
    % Input: face list and the current point
    % Output: indices of faces in the face list that are adjacent to the current point
    %         (can also be obtained using a point-face adjacency table)
    face_num=size(face,1);
    neigh_face_id=zeros(face_num,1);
    count=0;
    for i=1:face_num
        if any(ismember(face(i,:),point))
            count=count+1;
            neigh_face_id(count,1)=i;
        end
    end
    neigh_face_id=neigh_face_id(1:count,1);
end