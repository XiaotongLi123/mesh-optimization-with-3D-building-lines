function [normal_angle,region_edge_inside]=regionFaceNormalAngle(point,face,region_face_idx)
    % Compute the distribution of angles between adjacent faces in the influence region (in degrees)
    region_face=face(region_face_idx,:);
    normal=getFacesNormal(point,region_face);
    [~,region_edge]=getRegionPE(region_face);
    region_edge_inside=region_edge(region_edge(:,3)==1,1:2);
    inside_edge_num=size(region_edge_inside,1);
    normal_angle=zeros(inside_edge_num,1);
    for i=1:inside_edge_num
        neigh_face_pair=findEdgeNeighFace(region_face,region_edge_inside(i,:));
        normal_angle(i,1)=getVecAngle(normal(neigh_face_pair(1,1),:),normal(neigh_face_pair(2,1),:))*180/pi;
    end
end