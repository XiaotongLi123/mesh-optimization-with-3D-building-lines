function [length_sum]=getRegionLength(point,face,region_face_idx)
    region_face=face(region_face_idx,:);
    face_num=size(region_face,1);
    length=zeros(face_num,1);
    for i=1:face_num
        point_1=point(region_face(i,1),:);
        point_2=point(region_face(i,2),:);
        point_3=point(region_face(i,3),:);
        length(i,1)=point2pointDist(point_1,point_2)+point2pointDist(point_2,point_3)+point2pointDist(point_3,point_1);
    end
    length_sum=sum(length);
end