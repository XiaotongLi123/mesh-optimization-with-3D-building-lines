function [growth_face_idx]=regionGrowth(face,region_face_idx,PFneighbor,k)
    growth_face_idx=region_face_idx;
    for i=1:k
        growth_face=face(growth_face_idx,:);
        [growth_face_point,~]=getRegionPE(growth_face);
        growth_face_idx=getPointNeighFace(growth_face_point(:,1),PFneighbor);
    end
end