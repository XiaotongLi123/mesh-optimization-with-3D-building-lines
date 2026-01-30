function [repaired_region_face_idx,error]=holeFilling(region_face_idx,face,PFneighbor)
    % Hole filling
    % Step 1: Perform iterative region growing until no internal boundaries remain
    % Step 2: Extract the largest connected component and fill the other components
    % Initialization: if multiple connected components exist, select the largest one
    [region_point,~]=getRegionPE(face(region_face_idx,:));
    
    error=0;
    loop=0;  
    while loop<30
        % Region growing
        enlarge_face_idx=getPointNeighFace(region_point(:,1),PFneighbor);
        [region_point,region_edge]=getRegionPE(face(enlarge_face_idx,:));
        region_outline=region_edge(region_edge(:,3)==0,1:2);
        connect_outline=findConnectEdge(region_outline);
        connect_outline_count=size(connect_outline,1);
        if connect_outline_count==1
            break;
        end
        loop=loop+1;
    end
    if loop==30
        error=1;
    end
    
    % Identify and extract the largest connected component in the mesh    
    grow_face_idx=zeros(size(enlarge_face_idx,1),1);
    grow_face_num=0;
    for i=1:size(enlarge_face_idx,1)
        if ~any(ismember(region_face_idx,enlarge_face_idx(i,1)))
            grow_face_num=grow_face_num+1;
            grow_face_idx(grow_face_num,1)=enlarge_face_idx(i,1);
        end
    end

    [grow_connect_region]=findConnectRegion(grow_face_idx(1:grow_face_num,1),face);
    [~,max_grow_connect_id]=max(grow_connect_region(:,1));
    
    region_face_num=size(region_face_idx,1);
    region_face_idx=[region_face_idx;zeros(size(face,1),1)];

    for i=1:size(grow_connect_region,1)
        if i==max_grow_connect_id
            continue;
        end
        region_face_idx(region_face_num+1:region_face_num+grow_connect_region(i,1))=grow_connect_region(i,2:grow_connect_region(i,1)+1)';
        region_face_num=region_face_num+grow_connect_region(i,1);
    end
    repaired_region_face_idx=region_face_idx(1:region_face_num,1);
end