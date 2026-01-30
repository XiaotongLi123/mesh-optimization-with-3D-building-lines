function [new_face]=facepointUpdate(old_face,map)
    % Update face vertex indices based on the map
    face_num=size(old_face,1);
    new_face=zeros(face_num,3);
    for i=1:face_num
        for j=1:3                   
            new_face(i,j)=map(old_face(i,j),1);
        end
    end
end