function [normal]=getFacesNormal(point,face)
    face_num=size(face,1);
    normal=zeros(face_num,3);
    for i=1:face_num
        this_tri=point(face(i,:)',:);
        normal(i,:)=getPerFaceNormal(this_tri);
    end
end

function [face_normal]=getPerFaceNormal(triangle)
    vector=zeros(2,3);
    vector(1,:)=triangle(2,:)-triangle(1,:);
    vector(2,:)=triangle(3,:)-triangle(2,:);
    face_normal=crossProduct(vector);
    if norm(face_normal)~=0
        face_normal=face_normal/norm(face_normal);
    end
end

