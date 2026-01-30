function [face_edge]=getFaceEdge(face)
    face_edge=zeros(3,2);
    face_edge(1,:)=[face(1,2),face(1,3)];
    face_edge(2,:)=[face(1,3),face(1,1)];
    face_edge(3,:)=[face(1,1),face(1,2)];
end