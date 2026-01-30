function [new_region_face_idx,error]=regionSegment(point,face,region_face_idx,base_region_face_idx,line,constrained_edge,interval)
    pts_num=size(point,1);
    region_face=face(region_face_idx,:);
    [region_point,region_edge]=getRegionPE(region_face);
    
    region_point_num=size(region_point,1);
    region_edge_num=size(region_edge,1);
    constrained_edge_num=size(constrained_edge,1);

    error=0;

    if constrained_edge_num==0
        new_region_face_idx=region_face_idx;
        return;
    end

    % Construct the adjacency matrix and an undirected graph for region_point
    % Use unit weights to represent adjacency relationships
    connectMat=zeros(pts_num,pts_num); % Adjacency matrix of the influence region
    for i=1:region_edge_num
        connectMat(region_edge(i,1),region_edge(i,2))=point2pointDist(point(region_edge(i,1),:),point(region_edge(i,2),:));
        connectMat(region_edge(i,2),region_edge(i,1))=connectMat(region_edge(i,1),region_edge(i,2));
    end

    % Extract constrained edges within the influence region (excluding boundary edges) and connect them
    region_constrained_edge=zeros(constrained_edge_num,2);
    count=0;
    for i=1:constrained_edge_num
        in_region_edge=region_edge(region_edge(:,3)==1,1:2);
        [is_exist_line,~]=isExistLine(in_region_edge,constrained_edge(i,:));
        if is_exist_line
            count=count+1;
            region_constrained_edge(count,:)=constrained_edge(i,:);
        end
    end
    region_constrained_edge=region_constrained_edge(1:count,:);
    
    region_constrained_edge_num=size(region_constrained_edge,1);
    if region_constrained_edge_num==0
        new_region_face_idx=region_face_idx;
        return;
    end

    connect_region_constrained_edge_ori=findConnectEdge(region_constrained_edge);
    connect_constrained_edge_ori_num=size(connect_region_constrained_edge_ori,1);
    connect_region_constrained_edge_ori=[connect_region_constrained_edge_ori;zeros(100,size(connect_region_constrained_edge_ori,2))];

    % Eliminate edges that create self-loops (head-to-tail connections)
    for i=1:connect_constrained_edge_ori_num
        for j=2:connect_region_constrained_edge_ori(i,1)+1
            if any(ismember(connect_region_constrained_edge_ori(i,2:j-1),connect_region_constrained_edge_ori(i,j)))
                idx=find(connect_region_constrained_edge_ori(i,2:connect_region_constrained_edge_ori(i,1)+1)==connect_region_constrained_edge_ori(i,j));
                if idx(1,2)==connect_region_constrained_edge_ori(i,1)&&idx(1,1)~=1
                    connect_constrained_edge_ori_num=connect_constrained_edge_ori_num+1;
                    connect_region_constrained_edge_ori(connect_constrained_edge_ori_num,1)=idx(1,2)-idx(1,1)+1;
                    connect_region_constrained_edge_ori(connect_constrained_edge_ori_num,2:connect_region_constrained_edge_ori(connect_constrained_edge_ori_num,1)+1)=connect_region_constrained_edge_ori(i,idx(1,1)+1:idx(1,2)+1);
                    connect_region_constrained_edge_ori(i,1)=connect_region_constrained_edge_ori(i,1)-connect_region_constrained_edge_ori(connect_constrained_edge_ori_num,1)+1;
                    connect_region_constrained_edge_ori(i,idx(1,1)+2:idx(1,2)+1)=0;
                end
                if idx(1,1)==1&&idx(1,2)~=connect_region_constrained_edge_ori(i,1)
                    connect_constrained_edge_ori_num=connect_constrained_edge_ori_num+1;
                    connect_region_constrained_edge_ori(connect_constrained_edge_ori_num,1)=connect_region_constrained_edge_ori(i,1)-(idx(1,2)-idx(1,1)+1)+1;
                    connect_region_constrained_edge_ori(connect_constrained_edge_ori_num,2:connect_region_constrained_edge_ori(connect_constrained_edge_ori_num,1)+1)=connect_region_constrained_edge_ori(i,idx(1,2)+1:connect_region_constrained_edge_ori(i,1)+1);
                    connect_region_constrained_edge_ori(i,idx(1,2)+2:connect_region_constrained_edge_ori(i,1)+1)=0;
                    connect_region_constrained_edge_ori(i,1)=idx(1,2)-idx(1,1)+1;
                end
                break;
            end
        end
    end
    connect_region_constrained_edge=connect_region_constrained_edge_ori(1:connect_constrained_edge_ori_num,:);
    old_connect_constrained_edge_num=size(connect_region_constrained_edge,1);
    connect_constrained_edge_num=old_connect_constrained_edge_num;
    coMat=zeros(100,100);% In the correlation matrix and shortest path search, self-connections are not allowed
    
    % Process closed connected components (split the closed constrained line into two segments, 
    % but only extend one of the segments)
    for i=1:old_connect_constrained_edge_num
        start_point=connect_region_constrained_edge(i,2);
        
        end_point=connect_region_constrained_edge(i,1+connect_region_constrained_edge(i,1));
        if start_point~=end_point
            % If a connected component is not closed, check for possible connections to other connected components
            for j=1:old_connect_constrained_edge_num
                if j==i
                    continue;
                end
                if any(ismember(connect_region_constrained_edge(j,2:connect_region_constrained_edge(j,1)+1),start_point))||...
                        any(ismember(connect_region_constrained_edge(j,2:connect_region_constrained_edge(j,1)+1),end_point))
                    coMat(i,j)=1;
                    coMat(j,i)=1;
                end
            end
            continue;
        end
        length=round((connect_region_constrained_edge(i,1)-1)/2);
        new_connect_region_constrained_edge=zeros(1,size(connect_region_constrained_edge,2));
        new_connect_region_constrained_edge(1,1)=connect_region_constrained_edge(i,1)-length;
        new_connect_region_constrained_edge(1,2:new_connect_region_constrained_edge(1,1)+1)=connect_region_constrained_edge(i,2+length:1+connect_region_constrained_edge(i,1));
        connect_region_constrained_edge(i,1)=length+1;
        connect_region_constrained_edge(i,3+length:end)=0; % Note: include the truncation points twice here to ensure the constrained line segments remain connected
        connect_constrained_edge_num=connect_constrained_edge_num+1;
        connect_region_constrained_edge=[connect_region_constrained_edge;new_connect_region_constrained_edge]; % Connected component obtained after splitting a closed connected component
        coMat(i,connect_constrained_edge_num)=1;
        coMat(connect_constrained_edge_num,i)=1;
    end

    coMat=coMat(1:connect_constrained_edge_num,1:connect_constrained_edge_num);
    for i=1:connect_constrained_edge_num
        coMat(i,i)=1; % The constrained line is connected to itself
    end

    % Update coMat using a graph traversal algorithm
    coMat=coMatClustering(coMat);

    % Shortest path search: determine the shortest path from each endpoint to other constrained line segments or boundaries, then extend
    % Important points:
    % 1) Extend along the shortest path in each iteration, then update
    % 2) Update the adjacency matrix in each iteration
    process_sign=zeros(2*old_connect_constrained_edge_num,1);

    while sum(process_sign)~=2*old_connect_constrained_edge_num
        edge_shortest_path=zeros(2*old_connect_constrained_edge_num,region_point_num+1);
        edge_p_size=zeros(2*old_connect_constrained_edge_num,1);
        for i=1:2*old_connect_constrained_edge_num
            if process_sign(i,1)==1
                % Endpoints that have already been extended
                edge_shortest_path(i,1)=Inf;
                continue;
            end
            % Construct a new connectivity graph
            connectMat_temp=connectMat;
            if mod(i,2)==1
                this_end_point=connect_region_constrained_edge(ceil(i/2),2);
            else
                this_end_point=connect_region_constrained_edge(ceil(i/2),connect_region_constrained_edge(ceil(i/2),1)+1);
            end
            for j=1:connect_constrained_edge_num
                if coMat(ceil(i/2),j)==1
                    for k=2:connect_region_constrained_edge(j,1)+1
                        if connect_region_constrained_edge(j,k)~=this_end_point
                            connectMat_temp(connect_region_constrained_edge(j,k),:)=0;
                            connectMat_temp(:,connect_region_constrained_edge(j,k))=0;
                        end
                    end
                end
            end

            G_temp=graph(connectMat_temp);
            % Compute the shortest path for each endpoint
            point_set=region_point(region_point(:,2)==0,1);
            for j=1:connect_constrained_edge_num
                if coMat(ceil(i/2),j)==0
                    point_set=[point_set;connect_region_constrained_edge(j,2:connect_region_constrained_edge(j,1)+1)'];              
                end
            end
            point_set=unique(point_set);
            candidate_num=size(point_set,1);
            point_shortest_path=zeros(candidate_num,region_point_num+1);
            point_p_size=zeros(candidate_num,1);
            for j=1:candidate_num
                [p,d]=shortestpath(G_temp,this_end_point,point_set(j,1));
                point_shortest_path(j,1)=d;
                point_shortest_path(j,2:size(p,2)+1)=p;
                point_p_size(j,1)=size(p,2);
                if size(p,1)==0
                    point_shortest_path(j,1)=1000000;
                end
            end
            % Select the shortest path among all endpoints for extension, and perform the corresponding updates
            [~,min_id]=min(point_shortest_path(:,1));
            edge_shortest_path(i,:)=point_shortest_path(min_id,:);
            edge_p_size(i,1)=point_p_size(min_id,1);
        end    
        [~,all_min_id]=min(edge_shortest_path(:,1));
        extend_end_point=edge_shortest_path(all_min_id,edge_p_size(all_min_id,1)+1);
        % Extend the constrained edges
        if edge_p_size(all_min_id,1)~=0
            if mod(all_min_id,2)==0
                extend_edge=[connect_region_constrained_edge(ceil(all_min_id/2),2:connect_region_constrained_edge(ceil(all_min_id/2),1)+1),...
                    edge_shortest_path(all_min_id,3:edge_p_size(all_min_id,1)+1)];
                connect_region_constrained_edge(ceil(all_min_id/2),1)=connect_region_constrained_edge(ceil(all_min_id/2),1)+edge_p_size(all_min_id,1)-1;
                connect_region_constrained_edge(ceil(all_min_id/2),2:connect_region_constrained_edge(ceil(all_min_id/2),1)+1)=extend_edge;
            else
                extend_edge=[flip(edge_shortest_path(all_min_id,3:edge_p_size(all_min_id,1)+1)),...
                    connect_region_constrained_edge(ceil(all_min_id/2),2:connect_region_constrained_edge(ceil(all_min_id/2),1)+1)];
                connect_region_constrained_edge(ceil(all_min_id/2),1)=connect_region_constrained_edge(ceil(all_min_id/2),1)+edge_p_size(all_min_id,1)-1;
                connect_region_constrained_edge(ceil(all_min_id/2),2:connect_region_constrained_edge(ceil(all_min_id/2),1)+1)=extend_edge;
            end
        end
        % Process the markers
        process_sign(all_min_id,1)=1;

        % Associate markers
        for i=1:connect_constrained_edge_num
            if any(ismember(connect_region_constrained_edge(i,2:connect_region_constrained_edge(i,1)+1),extend_end_point))
                coMat(ceil(all_min_id/2),i)=1;
                coMat(i,ceil(all_min_id/2))=1;
            end
        end
        coMat=coMatClustering(coMat);
    end
    
    % Region clipping
    % Initialization
    extend_region_constrained_edge=zeros(region_edge_num,2);
    count=0;
    for i=1:connect_constrained_edge_num
        for j=2:connect_region_constrained_edge(i,1)
            count=count+1;
            extend_region_constrained_edge(count,:)=[connect_region_constrained_edge(i,j),connect_region_constrained_edge(i,j+1)];
        end
    end
    extend_region_constrained_edge=extend_region_constrained_edge(1:count,:);
    [region_segment,this_error]=segmentRegionGrowth(face,region_face_idx,extend_region_constrained_edge);
    
    if this_error==1
        error=1;
        new_region_face_idx=region_face_idx;
        return;
    end

    % Select the new LSM as the one with the shortest average distance 
    % from the clipped nodes to the constrained line segments
    region_num=size(region_segment,1);
    dist=zeros(region_num,1);
    for i=1:region_num
        this_region_idx=region_segment(i,2:region_segment(i,1)+1)';
        [~,~,dist(i,1),~]=line2meshDist(point,face,this_region_idx,line,interval);
    end
    [~,min_ind]=myMin(dist);
    if size(min_ind,1)==1
        new_region_face_idx=region_segment(min_ind,2:region_segment(min_ind,1)+1)';
    else
        area=zeros(size(min_ind,1),1);
        for i=1:size(min_ind,1)
            this_region_idx=region_segment(min_ind(i,1),2:region_segment(min_ind(i,1),1)+1)';
            area(i,1)=getRegionArea(point,face,this_region_idx);
        end
        [~,min_area_ind]=min(area);
        new_region_face_idx=region_segment(min_ind(min_area_ind,1),2:region_segment(min_ind(min_area_ind,1),1)+1)';
    end

    % LSM expansion to prevent incomplete regions after clipping
    new_region_face_idx=extendSegmentRegion(point,face,region_face_idx,base_region_face_idx,new_region_face_idx,line,region_constrained_edge,interval);
end


function [extend_region_segment]=extendSegmentRegion(point,face,region_face_idx,base_region_face_idx,segment_region_face_idx,line,constrained_edge,interval)
    [~,dist]=line2meshDist(point,face,segment_region_face_idx,line,interval);
    sign=1;
    while sign==1
        sign=0;
        region_face=face(segment_region_face_idx,:);
        [outline_point_sort,outline_sort]=getOutlineSort(region_face);

        % Extract constrained line segments on the boundary
        constrained_edge_num=size(constrained_edge,1);
        segment_region_constrained_edge=zeros(constrained_edge_num,2);
        count=0;
        for i=1:constrained_edge_num
            is_exist_line=isExistLine(outline_sort,constrained_edge(i,:));
            if is_exist_line
                count=count+1;
                segment_region_constrained_edge(count,:)=constrained_edge(i,:);
            end            
        end
        segment_region_constrained_edge=segment_region_constrained_edge(1:count,:);

        outline_num=size(outline_sort,1);
        
        % Extract candidate faces
        candidate_face_idx=zeros(outline_num,1);
        this_new_dist=zeros(outline_num,1);
        candidate_count=0;
        for i=1:outline_num
            is_exist_line=isExistLine(constrained_edge,outline_sort(i,:));
            % Cannot extend from the constrained edges
            if is_exist_line
                continue;
            end
            % Cannot extend if there are no adjacent faces
            neigh_face_idx=findEdgeNeighFace(face,outline_sort(i,:));
            this_candidate_face_idx=setdiff(intersect(neigh_face_idx,region_face_idx),segment_region_face_idx);
            if size(this_candidate_face_idx,1)==0||size(this_candidate_face_idx,2)==0
                continue;
            end
            % After extension, constrained edges must not be included again
            isExtend=1;
            candidate_face=face(this_candidate_face_idx,:);
            candidate_face_edge=getFaceEdge(candidate_face);
            for j=1:3
                is_exist_line=isExistLine(segment_region_constrained_edge,candidate_face_edge(j,:));
                if is_exist_line==1
                    isExtend=0;
                    break;
                end
            end
            if isExtend==0
                continue;
            end

            % After extension: 
            % If the distance remains unchanged, no holes are allowed to form; 
            % if the distance decreases, holes may be created
            if any(ismember(this_candidate_face_idx,base_region_face_idx))
                new_segment_region_face_idx=[segment_region_face_idx;this_candidate_face_idx];
                [~,temp_new_dist]=line2meshDist(point,face,new_segment_region_face_idx,line,interval);
                if temp_new_dist>dist
                    isExtend=0;
                end
                if temp_new_dist==dist
                    count=0;
                    for j=1:3
                        if any(ismember(outline_point_sort,candidate_face(1,j)))
                            count=count+1;
                        end
                    end
                    if count==3 
                        % If all three vertices of a new face lie on the boundary, 
                        % and two of its edges are not part of the original boundary, 
                        % a hole would be created; reject this case
                        is_exist_line_1=isExistLine(outline_sort,candidate_face_edge(1,:));
                        is_exist_line_2=isExistLine(outline_sort,candidate_face_edge(2,:));
                        is_exist_line_3=isExistLine(outline_sort,candidate_face_edge(3,:));
                        if is_exist_line_1+is_exist_line_2+is_exist_line_3==1
                            isExtend=0;
                        end
                    end
                end
                if isExtend==1
                    candidate_count=candidate_count+1;
                    candidate_face_idx(candidate_count,1)=this_candidate_face_idx;
                    this_new_dist(candidate_count,1)=temp_new_dist;
                end
            end
        end
        if candidate_count~=0
            sign=1;
            candidate_dist=zeros(candidate_count,1);
            for i=1:candidate_count
                tri=point(face(candidate_face_idx(i,1),:)',:);
                candidate_dist(i,1)=0.5*line2faceDist(tri,line,interval)+0.5*line2faceDist2(tri,line);
            end
            [~,add_face_id]=min(candidate_dist);
            add_face_idx=candidate_face_idx(add_face_id,1);
            dist=this_new_dist(add_face_id,1);
            segment_region_face_idx=[segment_region_face_idx;add_face_idx];
        end
    end
    extend_region_segment=segment_region_face_idx;
end

function [region_segment,error]=segmentRegionGrowth(face,region_face_idx,extend_region_constrained_edge)
    error=0;
    region_face=face(region_face_idx,:);
    region_face_num=size(region_face_idx,1);
    [~,region_edge]=getRegionPE(region_face);
    region_edge_num=size(region_edge,1);
    outline_edge=region_edge(region_edge(:,3)==0,1:2);
    outline_edge_num=size(outline_edge,1);

    % Initialization
    used_face=zeros(region_face_num,1);
    extend_region_constrained_edge_num=size(extend_region_constrained_edge,1);

    loop=0;
    region_segment=zeros(region_face_num,region_face_num+1);
    region_segment_count=0;
    while sum(used_face)~=region_face_num&&loop<10000
        this_segment_region_face_idx=zeros(region_face_num,1);
        % Update edge usage after each growth step
        used_edge=zeros(region_edge_num,1);
        for i=1:outline_edge_num  % Boundary edges are not allowed to grow
            [~,exist_line_id]=isExistLine(region_edge(:,1:2),outline_edge(i,:));
            used_edge(exist_line_id,1)=used_edge(exist_line_id,1)+1;
        end
        for i=1:extend_region_constrained_edge_num  % Extended constrained edges are not allowed to grow
            [~,exist_line_id]=isExistLine(region_edge(:,1:2),extend_region_constrained_edge(i,:));
            used_edge(exist_line_id,1)=used_edge(exist_line_id,1)+1;
        end

        % Seed triangle
        idx=(1:region_face_num)';
        candidate_face_idx=idx(used_face==0,1);
        seed_face_idx=candidate_face_idx(1,1);
        seed_edge=getFaceEdge(region_face(seed_face_idx,:));
        for i=1:3
            [~,exist_line_id]=isExistLine(region_edge(:,1:2),seed_edge(i,:));
            used_edge(exist_line_id,1)=used_edge(exist_line_id,1)+1;
        end
        
        this_segment_region_face_idx(1,1)=seed_face_idx;
        count=1;
        loop_start=count;
        loop_end=count;
        % growing
        while loop_start<=loop_end
            for i=loop_start:loop_end
                this_face=region_face(this_segment_region_face_idx(i,1),:);
                this_face_edge=getFaceEdge(this_face);
                for j=1:3
                    [~,exist_line_id]=isExistLine(region_edge(:,1:2),this_face_edge(j,:));
                    % If an edge has been used twice, it is no longer allowed to grow
                    if used_edge(exist_line_id,1)==2
                        continue;
                    end
                    % growing
                    third_face_id=findCertainNeighFace(region_face,(1:region_face_num)',this_face_edge(j,:),this_segment_region_face_idx(i,1));
                    if ~any(ismember(this_segment_region_face_idx(1:count,1),third_face_id))
                        count=count+1;
                        this_segment_region_face_idx(count,1)=third_face_id;
                        growth_edge=getFaceEdge(region_face(third_face_id,:));
                        for k=1:3
                            [~,exist_line_id]=isExistLine(region_edge(:,1:2),growth_edge(k,:));
                            used_edge(exist_line_id,1)=used_edge(exist_line_id,1)+1;
                        end
                    end
                end
                used_face(this_segment_region_face_idx(i,1),1)=1;
            end
            loop_start=loop_end+1;
            loop_end=count;
            loop=loop+1;
        end
        this_segment_region_face_idx=[count;region_face_idx(this_segment_region_face_idx(1:count),1);zeros(region_face_num-count,1)]';
        region_segment_count=region_segment_count+1;
        region_segment(region_segment_count,:)=this_segment_region_face_idx;
    end
    region_segment=region_segment(1:region_segment_count,:);
    if loop==10000
        region_segment=[];
        error=1;
    end
end

% Return all indices of the minimum elements in a column vector
function [min_val,min_id]=myMin(vec)
    min_val=min(vec);
    min_id=zeros(size(vec,1),1);
    count=0;
    for i=1:size(vec,1)
        if vec(i,1)<=1.005*min_val
            count=count+1;
            min_id(count,1)=i;
        end
    end
    min_id=min_id(1:count,1);
end

% Update coMat using a graph traversal algorithm
function [new_coMat]=coMatClustering(coMat) 
    % Aggregate adjacent line segments
    G=graph(coMat);

    line_num=size(coMat,1);
    line_sign=zeros(line_num,1);

    cluster_ori=zeros(line_num,line_num+1);
    cluster_num=0;
    
    for i=1:line_num      
        if line_sign(i,1)==1
            continue;
        end
        
        V=bfsearch(G,i);
        line_sign(V,1)=1;
        
        cluster_num=cluster_num+1;
        cluster_ori(cluster_num,1)=size(V,1);
        cluster_ori(cluster_num,2:cluster_ori(cluster_num,1)+1)=V';
    end

    new_coMat=zeros(line_num,line_num);
    for i=1:cluster_num
        for j=2:cluster_ori(i,1)+1
            new_coMat(cluster_ori(i,j),cluster_ori(i,2:cluster_ori(i,1)+1))=1;
        end
    end
end

