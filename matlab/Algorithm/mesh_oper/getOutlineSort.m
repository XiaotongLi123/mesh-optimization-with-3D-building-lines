function [outline_point_sort,outline_sort]=getOutlineSort(region_face)
    % Extract the boundary of LSM
    % Note: if faces are stored in counter-clockwise order (1,2 → 2,3 → 3,1),
    % the extracted boundary edges will also follow the counter-clockwise order
    [~,region_edge]=getRegionPE(region_face);

    outline=zeros(size(region_edge,1)-sum(region_edge(:,3)),2);
    outline_num=0;
    for i=1:size(region_edge,1)
        if region_edge(i,3)==0
            outline_num=outline_num+1;
            outline(outline_num,:)=region_edge(i,1:2);
        end
    end
   
    % Extract the boundary edges in sequential order
    outline_sort=regionEdgeSort(outline);
    outline_point_sort=outline_sort(:,1);
end

function [outline_sort]=regionEdgeSort(outline)
    outline_num=size(outline,1);
    for i=1:outline_num-1
        if outline(i,2)==outline(i+1,1)
            continue;
        end
        for j=i+2:outline_num
            if outline(i,2)==outline(j,1)
                t=outline(j,:);
                outline(j,:)=outline(i+1,:);
                outline(i+1,:)=t;
            end
        end
    end
    outline_sort=outline;
end