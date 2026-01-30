function [third_face_id]=findCertainNeighFace(face,face_idx,line,this_face_idx)
    % Find the face adjacent to a given edge of the current face
    % Input: faces in the influence region, the current face, and one of its edges
    % Output: ID of the face adjacent to the current edge; returns 0 if no such face exists
    third_face_id=zeros(1,1);
    for i=1:size(face,1)
        if any(ismember(face(i,:),line(1,1)))&&any(ismember(face(i,:),line(1,2)))
            if face_idx(i,1)~=this_face_idx
                third_face_id(1,1)=i;
                break;
            end
        end
    end
end