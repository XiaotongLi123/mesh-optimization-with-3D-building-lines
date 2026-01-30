function [max_dist,mean_dist,median_dist,node_dist]=line2meshDist(point,face,region_face_idx,line,interval)
    % Distance from a line to a mesh
    node=getNode(line,interval);
    node_num=size(node,1);

    p2m_dist=zeros(node_num,1);
    for i=1:node_num
        p2m_dist(i,1)=point2meshDist(point,face,region_face_idx,node(i,:));
    end
    node_dist=p2m_dist;
    max_dist=max(p2m_dist);
    mean_dist=mean(p2m_dist);
    median_dist=median(p2m_dist);
end
