function [l2f_dist]=line2faceDist(tri,line,interval)
    % Spatial distance from a line to a triangular face
    node=getNode(line,interval);
    node_num=size(node,1);

    p2f_dist=zeros(node_num,1);
    for i=1:node_num
        p2f_dist(i,1)=point2faceDist(node(i,:),tri);
    end
    l2f_dist=min(p2f_dist);
end