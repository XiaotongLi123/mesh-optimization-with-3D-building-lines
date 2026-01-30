function [new_region_face_idx,error]=regionConnect(point,face,region_face_idx,PFneighbor)
    % Multi-connected influence region connection algorithm
    % First, perform region growing until the influence regions are connected
    % Then, construct a connectivity graph for the grown regions (simplified graph to accelerate the algorithm),
    % use the shortest path algorithm to find the shortest path, and grow along this path
    ori_connect_region=findConnectRegion(region_face_idx,face);
    ori_connect_region_count=size(ori_connect_region,1);
    connect_region_count=ori_connect_region_count;
    error=0;
    if connect_region_count==1||connect_region_count==0
        new_region_face_idx=region_face_idx;
        return;
    end

    % ------------Region growth------------
    loop=0;
    enlarge_face_idx=region_face_idx;
    while connect_region_count~=1&&loop<30
        enlarge_region_face=face(enlarge_face_idx,:);
        [enlarge_region_point,~]=getRegionPE(enlarge_region_face);
        enlarge_face_idx=getPointNeighFace(enlarge_region_point(:,1),PFneighbor);
        connect_region=findConnectRegion(enlarge_face_idx,face);
        connect_region_count=size(connect_region,1);
        loop=loop+1;
    end
    if loop==30
        error=1;
        new_region_face_idx=region_face_idx;
        return;
    end

    % ------------Shortest-path based growth------------
    % Construct the adjacency matrix and undirected graph for region_point, using edge lengths as weights
    pts_num=size(point,1);
    enlarge_region_face=face(enlarge_face_idx,:);
    [~,enlarge_region_edge]=getRegionPE(enlarge_region_face);
    region_edge_num=size(enlarge_region_edge,1);
    connectMat=zeros(pts_num,pts_num);
    weight=zeros(region_edge_num,1);
    for i=1:region_edge_num
        weight(i,1)=point2pointDist(point(enlarge_region_edge(i,1),:),point(enlarge_region_edge(i,2),:));
        connectMat(enlarge_region_edge(i,1),enlarge_region_edge(i,2))=weight(i,1);
        connectMat(enlarge_region_edge(i,2),enlarge_region_edge(i,1))=weight(i,1);
    end
    G=graph(connectMat);
    
    % Shortest path search:
    % 1. Search for the shortest paths between all pairs of connected components
    % 2. In each iteration, connect the two components with the shortest path and update the connected components list
    loop=0;
    while ori_connect_region_count~=1&&loop<30
        % Extract the boundary points for each connected component
        outline_all=zeros(size(ori_connect_region,1),1000);
        for i=1:ori_connect_region_count
            this_region_face_idx=ori_connect_region(i,2:ori_connect_region(i,1)+1);
            this_region_point=getRegionPE(face(this_region_face_idx,:));
            outline=this_region_point(this_region_point(:,2)==0,1);
            outline_all(i,1)=size(outline,1);
            outline_all(i,2:outline_all(i,1)+1)=outline';
        end

        shortest_path=zeros(ori_connect_region_count*ori_connect_region_count,1+pts_num);
        shortest_path_connect_region=zeros(ori_connect_region_count*ori_connect_region_count,2);
        % Shortest path between regions
        for i=1:ori_connect_region_count
            for j=1:ori_connect_region_count
                if i>=j
                    continue;
                end
                % Shortest path between the boundary points of two regions
                this_shortest_path=zeros(outline_all(i,1)*outline_all(j,1),1+pts_num);
                d_all=zeros(outline_all(i,1)*outline_all(j,1),1);
                p_size=zeros(outline_all(i,1)*outline_all(j,1),1);
                for m=1:outline_all(i,1)
                    for n=1:outline_all(j,1)
                        [p,d_all((m-1)*outline_all(j,1)+n,1)]=shortestpath(G,outline_all(i,1+m),outline_all(j,1+n));
                        this_shortest_path((m-1)*outline_all(j,1)+n,1:size(p,2))=p;
                        p_size((m-1)*outline_all(j,1)+n,1)=size(p,2);
                    end
                end
                [~,min_id]=min(d_all);
                shortest_path((i-1)*ori_connect_region_count+j,1)=p_size(min_id);
                shortest_path((i-1)*ori_connect_region_count+j,2:p_size(min_id)+1)=this_shortest_path(min_id,1:p_size(min_id));
                shortest_path_connect_region((i-1)*ori_connect_region_count+j,1)=i;
                shortest_path_connect_region((i-1)*ori_connect_region_count+j,2)=j;
            end
        end
        for i=1:ori_connect_region_count*ori_connect_region_count
            if shortest_path(i,1)==0
                shortest_path(i,1)=Inf;
            end
        end
        [min_path_ptsnum,min_id]=min(shortest_path(:,1));
        neigh_face_idx=getPointNeighFace(shortest_path(min_id,2:min_path_ptsnum+1)',PFneighbor);
        start_region=shortest_path_connect_region(min_id,1);
        end_region=shortest_path_connect_region(min_id,2);
        new_face_idx=union(union(ori_connect_region(start_region,2:ori_connect_region(start_region,1)+1)',...
            ori_connect_region(end_region,2:ori_connect_region(end_region,1)+1)'),neigh_face_idx);
        ori_connect_region(start_region,1)=size(new_face_idx,1);
        ori_connect_region(start_region,2:size(new_face_idx,1)+1)=new_face_idx';
        ori_connect_region(end_region:end-1,:)=ori_connect_region(end_region+1:end,:);
        ori_connect_region_count=ori_connect_region_count-1;
        ori_connect_region=ori_connect_region(1:ori_connect_region_count,:);
    end
    if loop==30
        error=1;
        new_region_face_idx=region_face_idx;
        return;
    end
    new_region_face_idx=ori_connect_region(1,2:ori_connect_region(1,1)+1)';
end
