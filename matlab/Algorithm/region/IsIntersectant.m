function [region_face_idx,near_point]=IsIntersectant(point,face,knnsearch_face_idx,line,interval)
    % Sample nodes along a line segment and find the nearest points for each sampled node

    node=getNode(line,interval);
    node_num=size(node,1);
    
    % For each node, find the nearest faces and determine whether they intersect
    knn_face=face(knnsearch_face_idx,:);
    knn_face_num=size(knn_face,1);
    isRegionFace=zeros(knn_face_num,1);
    near_point=zeros(node_num,3);

    for n=1:node_num
        [~,nearest_face_idx]=point2meshDist(point,face,knnsearch_face_idx,node(n,:));
        
        nearest_face_num=size(nearest_face_idx,1);
        for i=1:nearest_face_num
            this_face_idx=knnsearch_face_idx(nearest_face_idx(i,1),1);  

            triangle=point(face(this_face_idx,:)',:);
            [~,nearest_point]=point2faceDist(node(n,:),triangle);

            near_point(n,:)=nearest_point;

            is_near=faceNearestPointPos(nearest_point,triangle);
            near_edge_num=sum(is_near);
            if near_edge_num==0
                isRegionFace(nearest_face_idx(i,1),1)=1;
                continue;
            end
            if near_edge_num==1
                neigh_face_id=findEdgeNeighFace(knn_face,face(this_face_idx,is_near==0));
                isRegionFace(neigh_face_id,1)=1;
                continue;
            end
            if near_edge_num==2
                neigh_face_id=findPointNeighFace(knn_face,face(this_face_idx,is_near==0));
                isRegionFace(neigh_face_id,1)=1;
                continue;
            end
        end
    end
    region_face_idx=knnsearch_face_idx(isRegionFace==1,1);
end

function [is_near]=faceNearestPointPos(proj_point,tri)
    is_near=zeros(3,1);
    dist_threshold=zeros(3,1);
    line=zeros(6,3);
    line(1:2,:)=tri(2:3,:);
    line(3:4,:)=[tri(1,:);tri(3,:)];
    line(5:6,:)=tri(1:2,:);
    for i=1:3
        dist_threshold(i,1)=0.25*point2lineDist(tri(i,:),line(2*i-1:2*i,:));
        if point2lineDist(proj_point,line(2*i-1:2*i,:))<=dist_threshold(i,1)
            is_near(i,1)=1;
        end
    end
end