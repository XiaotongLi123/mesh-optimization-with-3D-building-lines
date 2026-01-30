function [neigh_face_id]=findEdgeNeighFace(face,line)
    % Find the triangular faces adjacent to the current line segment
    % Input: face list and the current line segment
    % Output: indices of faces in the face list that are adjacent to the current line segment
    neigh_face_id=zeros(2,1);
    count=0;
    for i=1:size(face,1)
        if any(ismember(face(i,:),line(1,1)))&&any(ismember(face(i,:),line(1,2)))
            count=count+1;
            neigh_face_id(count,1)=i;
        end
        if count==2
            break;
        end
    end
    neigh_face_id=neigh_face_id(1:count,1);
end