function [have]=haveConstrainedEdge(face,region_face_idx,constrained_edge)
    region_face=face(region_face_idx,:);
    [~,region_edge]=getRegionPE(region_face);
    
    constrained_edge_num=size(constrained_edge,1);

    have=0;
    if constrained_edge_num==0
        return;
    end

    % Identify constrained edges within the influence region, excluding boundary edges, and link them
    for i=1:constrained_edge_num
        in_region_edge=region_edge(region_edge(:,3)==1,1:2);
        [is_exist_line,~]=isExistLine(in_region_edge,constrained_edge(i,:));
        if is_exist_line
            have=1;
        end
    end   
end