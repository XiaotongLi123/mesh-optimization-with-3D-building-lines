function [area]=getRegionArea(point,face,region_face_idx)
    face_num=size(region_face_idx,1);
    area=0;
    for i=1:face_num
        tri=point(face(region_face_idx(i,1),:)',:);
        area=area+getTriArea(tri);
    end
end