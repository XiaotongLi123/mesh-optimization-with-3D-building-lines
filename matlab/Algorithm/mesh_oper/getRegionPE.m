function [region_point,region_edge]=getRegionPE(region_face)
    % Extract points and edges within a region, and label them
    % Boundary points and edges are labeled as 0, internal points and edges are labeled as 1
    % For non-manifold structures, edge labels are typically greater than 1
    % region_point is a 2-column matrix, region_edge is a 3-column matrix
    region_face_num=size(region_face,1);
    region_point_ori=[zeros(3*region_face_num,1),ones(3*region_face_num,1)];
    region_edge_ori=zeros(3*region_face_num,3);
    
    region_edge_num=0;
    region_point_num=0;

    % Extract edges within the region
    for i=1:region_face_num
        this_face=region_face(i,:);
        for j=1:3
            if j==1
                this_edge=[this_face(1,1),this_face(1,2)];
            end
            if j==2
                this_edge=[this_face(1,2),this_face(1,3)]; 
            end
            if j==3
                this_edge=[this_face(1,3),this_face(1,1)];
            end
            [is_exist_edge,exist_line_id]=isExistLine(region_edge_ori(1:region_edge_num,1:2),this_edge);
            if is_exist_edge==1
                region_edge_ori(exist_line_id,3)=region_edge_ori(exist_line_id,3)+1;
            end
            if is_exist_edge==0
                region_edge_num=region_edge_num+1;
                region_edge_ori(region_edge_num,1:2)=this_edge;
            end
        end
    end
    region_edge=region_edge_ori(1:region_edge_num,:);

    % Extract vertices within the region
    for i=1:region_edge_num
        for j=1:2
            if ~any(ismember(region_point_ori(1:region_point_num,1),region_edge(i,j)))
                region_point_num=region_point_num+1;
                region_point_ori(region_point_num,1)=region_edge(i,j);
                if region_edge(i,3)==0
                    region_point_ori(region_point_num,2)=0;
                end
            else
                if region_edge(i,3)==0
                    region_point_ori(region_point_ori(:,1)==region_edge(i,j),2)=0;
                end
            end
        end
    end
    region_point=region_point_ori(1:region_point_num,:);
end