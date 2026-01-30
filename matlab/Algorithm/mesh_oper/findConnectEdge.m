function [connect_edge]=findConnectEdge(edge)
    % Find the connected components of connected edges
    % If the start and end points of a component are the same, the component is closed
    % Handle cases where a single point connects to multiple edges (currently consider up to three connections per point)
    edge_num=size(edge,1);
    max_length=edge_num+1;
    used=zeros(edge_num,1);
    connect_edge_count=0;

    connect_edge_ori=zeros(edge_num,max_length+1);
    for i=1:edge_num
        if used(i,1)==1
            continue;
        end
        used(i,1)=1;
        
        this_connect=zeros(1,max_length);
        left_neigh=-1;
        this_connect(1,end)=edge(i,1);
        this_connect(1,1)=edge(i,2);
        left_neigh_num=0;

        while left_neigh~=0
            left_neigh=0;
            for j=1:edge_num
                if used(j,1)==1
                    continue;
                end
                for k=1:2
                    if edge(j,k)==this_connect(1,end-left_neigh_num)
                        left_neigh_num=left_neigh_num+1;
                        this_connect(1,end-left_neigh_num)=edge(j,3-k);
                        used(j,1)=1;
                        left_neigh=1;
                        break;
                    end
                end
            end
        end
        
        right_neigh=-1;
        right_neigh_num=0;
        while right_neigh~=0
            right_neigh=0;
            for j=1:edge_num
                if used(j,1)==1
                    continue;
                end
                for k=1:2
                    if edge(j,k)==this_connect(1,1+right_neigh_num)
                        right_neigh_num=right_neigh_num+1;
                        this_connect(1,1+right_neigh_num)=edge(j,3-k);
                        used(j,1)=1;
                        right_neigh=1;
                        break;
                    end
                end
            end
        end
        connect_edge_count=connect_edge_count+1;
        connect_edge_ori(connect_edge_count,1)=2+left_neigh_num+right_neigh_num;
        connect_edge_ori(connect_edge_count,2:left_neigh_num+2)=this_connect(1,end-left_neigh_num:end);
        connect_edge_ori(connect_edge_count,left_neigh_num+3:connect_edge_ori(connect_edge_count,1)+1)...
            =this_connect(1,1:1+right_neigh_num);
    end
    connect_edge=connect_edge_ori(1:connect_edge_count,:);
end