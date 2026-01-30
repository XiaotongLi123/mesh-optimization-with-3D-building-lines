function [p2m_dist,nearest_face_idx]=point2meshDist(point,face,region_face_idx,this_point)
    % Compute the distance from a point to a mesh, and return the closest triangular face
    region_face_num=size(region_face_idx,1);
    p2f_dist=zeros(region_face_num,1);
    for i=1:region_face_num
        p2f_dist(i,1)=point2faceDist(this_point,point(face(region_face_idx(i,1),:)',:));
    end
    [p2m_dist]=min(p2f_dist);
    nearest_face_idx=zeros(region_face_num,1);
    count=0;
    for i=1:region_face_num
        if p2f_dist(i,1)==p2m_dist
            count=count+1;
            nearest_face_idx(count,1)=i;
        end
    end
    nearest_face_idx=nearest_face_idx(1:count,1);
end