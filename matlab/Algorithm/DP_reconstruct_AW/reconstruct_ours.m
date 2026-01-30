function [new_face,new_point,normal_angle,constrained_line,opt_face_c,error]=reconstruct_ours(input_face,input_point,input_region_face_idx,line,PFneighbor,constrained_line,interval)   
    % Mesh reconstruction with 3D line-element constraints
    % Input: vertices and faces of the full mesh, indices of faces in the influence region,
    %        constrained line elements
    % Output: reconstructed vertices and faces of the full mesh,
    %         Hausdorff distance,
    %         dihedral angles within the influence region (for experimental analysis),
    %         constrained line segments (preserved during subsequent reconstruction),
    %         reconstruction success flag
    
    input_face_num=size(input_face,1); % Number of input faces
    input_region_face_num=size(input_region_face_idx,1); % Number of faces in LSM

    new_face_ori=zeros(input_face_num+10*input_region_face_num,3);
    new_face_ori(1:input_face_num,:)=input_face; % Construct a new face list
    error=0;

    %---------------Extract influence-region boundaries and compute reconstruction weights---------------
    % Extract the boundary of the influence region in order
    % Note: if faces are stored in counterclockwise order (1,2–2,3–3,1),
    % the extracted boundary is also counterclockwise. Therefore, boundary
    % consistency can be checked at the end to determine whether a normal
    % flip is required for the reconstructed mesh.
    region_face=input_face(input_region_face_idx,:); % Region before reconstruction
    [region_point,~]=getRegionPE(region_face); % Extract all vertices in the influence region
    [edge_point_sort,edge_sort]=getOutlineSort(region_face); % Extracted boundary vertices arranged in counterclockwise order
    outline_num=size(edge_point_sort,1);
    
    % If a constrained line fully coincides with the boundary,
    % use the boundary directly as the constrained line element
    [constrained_outline]=findConstrainedOutline(input_point,edge_point_sort,line);
    if size(constrained_outline,1)~=0
        new_face=input_face;
        new_point=input_point;
        constrained_line=[constrained_line;constrained_outline];
        normal_angle=zeros(1,3);
        return;
    end

    enlarge_face_idx=getPointNeighFace(region_point(:,1),PFneighbor); % Find the triangle adjacent to each boundary edge and compute its normal vector
    neigh_face_normal=zeros(outline_num,3);
    is_constrained=zeros(outline_num,1);
    for i=1:outline_num
        edge_neigh_face_id=findEdgeNeighFace(input_face(enlarge_face_idx,:),edge_sort(i,:));
        edge_neigh_face_idx=enlarge_face_idx(edge_neigh_face_id,:);
        if size(edge_neigh_face_idx,1)==1
            continue;
        end
        neigh_face_idx=setdiff(edge_neigh_face_idx,input_region_face_idx);
        neigh_face=input_face(neigh_face_idx,:);
        is_constrained(i,1)=isExistLine(constrained_line,edge_sort(i,:));
        neigh_face_normal(i,:)=getFacesNormal(input_point,neigh_face);
    end
    %---------------Split constrained line elements---------------
    node=getNode(line,interval); 
    [weight_1,weight_2]=getWeight(input_point,input_face,input_region_face_idx,0.5,0.5); % Weight computation function
    weight=weight_2/weight_1;

    %---------------Mesh retriangulation using dynamic programming---------------
    % Reconstruct the mesh using a dynamic programming algorithm
    % When incorporating constrained line elements, the endpoints of the constraints
    % are connected to two boundary vertices each, and exhaustive search is used to find the optimal connection
    % Once the connected endpoints are determined, the reconstruction problem
    % is reduced to two 3D polygon reconstruction problems, which are solved using dynamic programming
    [opt_face,opt_val]=DP_reconstruct_ours(input_point(edge_point_sort,:),node,neigh_face_normal,is_constrained,weight);
    if opt_val==inf % Retriangulation failed
        error=1;
        new_face=input_face;
        new_point=input_point;
        normal_angle=zeros(1,3);
        return;
    end

    %---------------Post-processing---------------
    % Construct a map from indices of points in the reconstructed region to indices in the full mesh,
    % and adjust face vertex indices accordingly for various cases
    % Add line segment endpoints and nodes to the point list; their indices range from (point_num - node_num) to point_num
    % Construct the list of constrained lines
    [map,new_constrained_line,new_point]=postProcess(input_point,node,edge_point_sort);
    constrained_line=[constrained_line;new_constrained_line];
    opt_face_c=faceChange(opt_face,map);

    reconstruct_face_num=size(opt_face_c,1);

    % Remove deleted triangles
    new_face_ori(input_region_face_idx,:)=0;
    new_face_num=input_face_num-input_region_face_num+reconstruct_face_num;
    new_face=zeros(new_face_num,3);
    count=0;
    for i=1:input_face_num
        if new_face_ori(i,1)~=0
            count=count+1;
            new_face(count,:)=new_face_ori(i,:);
        end
    end
    new_face(count+1:new_face_num,:)=opt_face_c;


    % Set all points inside the LSM to zero without removing them from the list,
    % keeping their original indices unchanged
    for i=1:size(region_point,1)
        if region_point(i,2)==1
            new_point(region_point(i,1),:)=0;
        end
    end

    %---------------Evaluation---------------
    % Analyze dihedral angles of faces in the influence region
    normal_angle=regionFaceNormalAngle(new_point,new_face,(new_face_num-reconstruct_face_num+1:new_face_num)');
end

function [new_face]=faceChange(old_face,map)
    % Update face vertex indices based on the map
    face_num=size(old_face,1);
    new_face=zeros(face_num,3);
    for i=1:face_num
        for j=1:3
            new_face(i,j)=map(old_face(i,j),1);
        end
    end
end

function [map,constrained_line,new_point]=postProcess(input_point,node,edge_point_sort)
    % Create a map, set of constrained line segments, and the new point list of the reconstructed mesh
    [is_exist_point_1,exist_point_id_1]=isExistPoint(input_point(edge_point_sort,:),node(1,:));
    [is_exist_point_2,exist_point_id_2]=isExistPoint(input_point(edge_point_sort,:),node(end,:));
    node_num=size(node,1);
    input_point_num=size(input_point,1);

    constrained_line=zeros(node_num-1,2);

    if is_exist_point_1==0
        if is_exist_point_2==0
            map=[edge_point_sort;(input_point_num+1:input_point_num+node_num)'];
            new_point=[input_point;node];

            for i=1:node_num-1
                constrained_line(i,:)=[input_point_num+i,input_point_num+i+1];
            end
        else
            map=[edge_point_sort;(input_point_num+1:input_point_num+node_num-1)'];
            new_point=[input_point;node(1:end-1,:)];

            for i=1:node_num-2
                constrained_line(i,:)=[input_point_num+i,input_point_num+i+1];
            end
            constrained_line(node_num-1,:)=[input_point_num+node_num-1,edge_point_sort(exist_point_id_2,1)];
        end
    else
        if is_exist_point_2==0
            map=[edge_point_sort;(input_point_num+1:input_point_num+node_num-1)'];
            new_point=[input_point;node(2:end,:)];

            constrained_line(1,:)=[edge_point_sort(exist_point_id_1,1),input_point_num+1];
            for i=1:node_num-2
                constrained_line(i+1,:)=[input_point_num+i,input_point_num+i+1];
            end
        else
            map=[edge_point_sort;(input_point_num+1:input_point_num+node_num-2)'];
            new_point=[input_point;node(2:end-1,:)];
            if node_num>2
                constrained_line(1,:)=[edge_point_sort(exist_point_id_1,1),input_point_num+1];
                for i=1:node_num-3
                    constrained_line(i+1,:)=[input_point_num+i,input_point_num+i+1];
                end
                constrained_line(node_num-1,:)=[input_point_num+node_num-2,edge_point_sort(exist_point_id_2,1)];
            else
                constrained_line(1,:)=[edge_point_sort(exist_point_id_1,1),edge_point_sort(exist_point_id_2,1)];
            end
        end
    end
end

function [constrained_outline]=findConstrainedOutline(point,outline_point_sort,line)
    outline_point_num=size(outline_point_sort,1);
    constrained_outline=[];
    left=0;
    right=0;
    for i=1:outline_point_num
        if point2pointDist(point(outline_point_sort(i,1),:),line(1,:))<1e-4
            left=i;
        end
        if point2pointDist(point(outline_point_sort(i,1),:),line(2,:))<1e-4
            right=i;
        end
    end
    overlapping_sign=ones(2,1);
    if left~=right&&left~=0&&right~=0
        for i=1:outline_point_num
            if i>=min(left,right)&&i<=max(left,right)
                if point2segmentDist(point(outline_point_sort(i,1),:),line)>1e-4
                    overlapping_sign(1,1)=0;
                end
            else
                if point2segmentDist(point(outline_point_sort(i,1),:),line)>1e-4
                    overlapping_sign(2,1)=0;
                end
            end
        end
    end
    if overlapping_sign(1,1)==1&&overlapping_sign(2,1)==0
        constrained_outline=zeros(max(left,right)-min(left,right),2);
        for i=min(left,right):max(left,right)-1
            constrained_outline(i+1-min(left,right),:)=[outline_point_sort(i,1),outline_point_sort(i+1,1)];
        end
    end
    if overlapping_sign(2,1)==1&&overlapping_sign(1,1)==0
        constrained_outline=zeros(outline_point_num-(max(left,right)-min(left,right)),2);
        for i=1:min(left,right)-1
            constrained_outline(i,:)=[outline_point_sort(i,1),outline_point_sort(i+1,1)];
        end
        for i=max(left,right):outline_point_num-1
            constrained_outline(i+min(left,right)-1,:)=[outline_point_sort(i,1),outline_point_sort(i+1,1)];
        end
        constrained_outline(outline_point_num-(max(left,right)-min(left,right)),:)=[outline_point_sort(outline_point_num,1),outline_point_sort(1,1)];
    end
end
