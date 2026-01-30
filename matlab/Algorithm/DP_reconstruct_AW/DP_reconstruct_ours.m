function [opt_face,opt_val]=DP_reconstruct_ours(point,node,neigh_face_normal,is_constrained,weight)
    % Influence-region reconstruction using dynamic programming,
    % where the polygonal boundary consists of multiple boundary points and a constrained line,
    % aiming to minimize the maximum dihedral angle of triangles through global optimization    

    %----------------Initialize----------------
    % Initialize the optimal value
    n=size(point,1);
    opt_face=[];
    opt_val=inf;
    %----------------Divide the influence region into two 3D polygons and reconstruct each using dynamic programming----------------
    % Handle four cases depending on whether the endpoints of the constrained line lie on the boundary of the influence region.
    % Endpoints not on the boundary are connected via exhaustive search to determine the optimal partition into two spatial polygons.
    % The region is then divided accordingly and reconstructed using dynamic programming.
    [is_exist_point_1,exist_point_id_1]=isExistPoint(point,node(1,:));
    [is_exist_point_2,exist_point_id_2]=isExistPoint(point,node(end,:));
    % Both endpoints lie on the boundary
    if is_exist_point_1==1&&is_exist_point_2==1
        left=exist_point_id_1(1,1);
        right=exist_point_id_2(1,1);
        node=node(2:end-1,:);
        % Compare the current reconstruction with the previous optimum and update if improved.
        % DP_polygonReconstruct is called twice within the update function to reconstruct the two spatial polygons.
        [opt_face,opt_val]=update(point,node,neigh_face_normal,is_constrained,left,right,opt_face,opt_val,weight,1,1);
        return;
    end

    KDmodel=KDTreeSearcher(point);% Use k-d tree to find the k nearest points to the left/right endpoints, reducing the number of iterations
    neigh=min(2,n);
    % Left endpoint is on the boundary; iterate over boundary points connected to the right endpoint
    if is_exist_point_1==1&&is_exist_point_2==0
        left=exist_point_id_1(1,1);
        index=knnsearch(KDmodel,node(end,:),'K',neigh);
        node=node(2:end,:);
        for i=1:neigh
            right=index(1,i);
            if right==left
                continue;
            end
            [opt_face,opt_val]=update(point,node,neigh_face_normal,is_constrained,left,right,opt_face,opt_val,weight,1,0);
            if opt_val~=Inf
                break;
            end
        end
        return;
    end
    % Right endpoint is on the boundary; iterate over boundary points connected to the left endpoint
    if is_exist_point_2==1&&is_exist_point_1==0
        right=exist_point_id_2(1,1);
        index=knnsearch(KDmodel,node(1,:),'K',neigh);
        node=node(1:end-1,:);

        for i=1:neigh
            left=index(1,i);
            if left==right
                continue;
            end
            [opt_face,opt_val]=update(point,node,neigh_face_normal,is_constrained,left,right,opt_face,opt_val,weight,0,1);
            if opt_val~=Inf
                break;
            end
        end
        return;
    end 
    % If neither endpoint lies on the boundary, iterate over boundary points connected to both endpoints
    index_1=knnsearch(KDmodel,node(1,:),'K',neigh);
    index_2=knnsearch(KDmodel,node(end,:),'K',neigh);
    for i=1:neigh
        for j=1:neigh
            left=index_1(1,i);
            right=index_2(1,j);
            if left==right
                continue;
            end
            [opt_face,opt_val]=update(point,node,neigh_face_normal,is_constrained,left,right,opt_face,opt_val,weight,0,0);
            if opt_val~=Inf
                break;
            end
        end
    end
end

function [opt_face,opt_val]=update(point,node,neigh_face_normal,is_constrained,left,right,opt_face,opt_val,weight,left_sign,right_sign)
    % Inputs: boundary vertices, nodes, left-endpoint connected boundary points, right-endpoint connected boundary points,
    %         previous reconstruction (triangles) and optimal value, initial normals, point groupings, two weights
    % Outputs: updated triangles and updated optimal value
    pts_num=size(point,1);
    node_num=size(node,1);
    n=pts_num;

    %-----------------Determine the points and map for each of the two spatial polygons-----------------
    % Since the boundary is stored in order (a linked list could be used in C++),
    % different cases need to be considered
    if left<right
        % Group the points
        point_1=[point(left:right,:);lineDirectionChange(node)];
        neigh_face_normal_1=zeros(size(point_1,1),3);
        neigh_face_normal_1(1:right-left,:)=neigh_face_normal(left:right-1,:);
        is_constrained_1=zeros(size(point_1,1),1);
        is_constrained_1(1:right-left,:)=is_constrained(left:right-1,:);
        point_2=[point(right:pts_num,:);point(1:left,:);node];
        neigh_face_normal_2=zeros(size(point_2,1),3);
        neigh_face_normal_2(1:pts_num-right+left,:)=[neigh_face_normal(right:pts_num,:);neigh_face_normal(1:left-1,:)];
        is_constrained_2=zeros(size(point_2,1),1);
        is_constrained_2(1:pts_num-right+left,:)=[is_constrained(right:pts_num,:);is_constrained(1:left-1,:)];
        % Establish correspondence between the new point list and the original point list
        map_1=[(left:right)';(n+node_num:-1:n+1)'];
        map_2=[(right:pts_num)';(1:left)';(n+1:n+node_num)'];
    else
        % Group the points
        point_1=[point(right:left,:);node];
        neigh_face_normal_1=zeros(size(point_1,1),3);
        neigh_face_normal_1(1:left-right,:)=neigh_face_normal(right:left-1,:);
        is_constrained_1=zeros(size(point_1,1),1);
        is_constrained_1(1:left-right,:)=is_constrained(right:left-1,:);
        point_2=[point(left:pts_num,:);point(1:right,:);lineDirectionChange(node)];
        neigh_face_normal_2=zeros(size(point_2,1),3);
        neigh_face_normal_2(1:pts_num-left+right,:)=[neigh_face_normal(left:pts_num,:);neigh_face_normal(1:right-1,:)];
        is_constrained_2=zeros(size(point_2,1),1);
        is_constrained_2(1:pts_num-left+right,:)=[is_constrained(left:pts_num,:);is_constrained(1:right-1,:)];
        % Establish correspondence between the new point list and the original point list
        map_1=[(right:left)';(n+1:n+node_num)'];
        map_2=[(left:pts_num)';(1:right)';(n+node_num:-1:n+1)'];
    end

 
    %-----------------Split into two polygonal optimization problems-----------------
    new_neigh_face_normal_1=neigh_face_normal_1;
    new_is_constrained_1=is_constrained_1;
    new_neigh_face_normal_2=neigh_face_normal_2;
    new_is_constrained_2=is_constrained_2;
    
    % First reconstruct polygon 1, then reconstruct polygon 2
    [opt_face_1,opt_val_1]=DP_PR_ours(point_1,neigh_face_normal_1,is_constrained_1,1,weight);
    opt_face_1_c=facepointUpdate(opt_face_1,map_1);

    if left<right
        for i=pts_num-right+left+1:size(point_2,1)
            if i~=size(point_2,1)
                edge_neigh_face_id=findEdgeNeighFace(opt_face_1_c,[map_2(i,1),map_2(i+1,1)]);
                if i~=pts_num-right+left+1
                    new_is_constrained_2(i,1)=1;
                else
                    if left_sign==1
                        new_is_constrained_2(i,1)=1;
                    end
                end
            else
                edge_neigh_face_id=findEdgeNeighFace(opt_face_1_c,[map_2(i,1),map_2(1,1)]);
                if right_sign==1
                    new_is_constrained_2(i,1)=1;
                end
            end    
            if size(edge_neigh_face_id,1)==0
                continue;
            end
            new_neigh_face_normal_2(i,:)=getFacesNormal([point;node],opt_face_1_c(edge_neigh_face_id,:));
        end
    else
        for i=pts_num-left+right+1:size(point_2,1)
            if i~=size(point_2,1)
                edge_neigh_face_id=findEdgeNeighFace(opt_face_1_c,[map_2(i,1),map_2(i+1,1)]);
                if i~=pts_num-left+right+1
                    new_is_constrained_2(i,1)=1;
                else
                    if right_sign==1
                        new_is_constrained_2(i,1)=1;
                    end
                end
            else
                edge_neigh_face_id=findEdgeNeighFace(opt_face_1_c,[map_2(i,1),map_2(1,1)]);
                if left_sign==1
                    new_is_constrained_2(i,1)=1;
                end
            end    
            if size(edge_neigh_face_id,1)==0
                continue;
            end
            new_neigh_face_normal_2(i,:)=getFacesNormal([point;node],opt_face_1_c(edge_neigh_face_id,:));
        end
    end

    [opt_face_2,opt_val_2]=DP_PR_ours(point_2,new_neigh_face_normal_2,new_is_constrained_2,1,weight);
    opt_face_2_c=facepointUpdate(opt_face_2,map_2);

    this_opt_val=opt_val_1+opt_val_2;
    opt_face_c=[opt_face_1_c;opt_face_2_c];
    normal_angle=regionFaceNormalAngle([point;node],opt_face_c,(1:size(opt_face_c,1))');

    if max(normal_angle)>75
        % Reconstruct polygon 2 first, then polygon 1 (if the first reconstruction contains folding)
        [opt_face_2,opt_val_2]=DP_PR_ours(point_2,neigh_face_normal_2,is_constrained_2,1,weight);
        opt_face_2_c=facepointUpdate(opt_face_2,map_2);
        
        if left<right
            for i=right-left+1:size(point_1,1)
                if i~=size(point_1,1)
                    edge_neigh_face_id=findEdgeNeighFace(opt_face_2_c,[map_1(i,1),map_1(i+1,1)]);
                    if i~=right-left+1
                        new_is_constrained_1(i,1)=1;
                    else
                        if right_sign==1
                            new_is_constrained_1(i,1)=1;
                        end
                    end
                else
                    edge_neigh_face_id=findEdgeNeighFace(opt_face_2_c,[map_1(i,1),map_1(1,1)]);
                    if left_sign==1
                        new_is_constrained_1(i,1)=1;
                    end
                end    
                if size(edge_neigh_face_id,1)==0
                    continue;
                end
                new_neigh_face_normal_1(i,:)=getFacesNormal([point;node],opt_face_2_c(edge_neigh_face_id,:));
            end
        else
            for i=left-right+1:size(point_1,1)
                if i~=size(point_1,1)
                    edge_neigh_face_id=findEdgeNeighFace(opt_face_2_c,[map_1(i,1),map_1(i+1,1)]);
                    if i~=left-right+1
                        new_is_constrained_1(i,1)=1;
                    else
                        if left_sign==1
                            new_is_constrained_1(i,1)=1;
                        end
                    end
                else
                    edge_neigh_face_id=findEdgeNeighFace(opt_face_2_c,[map_1(i,1),map_1(1,1)]);
                    if right_sign==1
                        new_is_constrained_1(i,1)=1;
                    end
                end    
                if size(edge_neigh_face_id,1)==0
                    continue;
                end
                new_neigh_face_normal_1(i,:)=getFacesNormal([point;node],opt_face_2_c(edge_neigh_face_id,:));
            end
        end
        [opt_face_1,opt_val_1]=DP_PR_ours(point_1,new_neigh_face_normal_1,new_is_constrained_1,1,weight);
        opt_face_1_c=facepointUpdate(opt_face_1,map_1);
    
        new_this_opt_val=opt_val_1+opt_val_2;
        if new_this_opt_val<this_opt_val
            opt_face_c=[opt_face_1_c;opt_face_2_c];
            this_opt_val=new_this_opt_val;
        end
    end

    %-----------------If the update condition is met, update the previously reconstructed result-----------------
    if this_opt_val<opt_val
        opt_val=this_opt_val;
        opt_face=opt_face_c;
    end
end

function [node_opposite]=lineDirectionChange(node)
    % Reverse the constrained line
    node_num=size(node,1);
    node_opposite=zeros(node_num,3);
    for i=1:node_num
        node_opposite(i,:)=node(node_num+1-i,:);
    end
end