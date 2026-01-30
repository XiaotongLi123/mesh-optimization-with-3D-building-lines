function [connect_region]=findConnectRegion(region_face_idx,face)
    % Divide the LSM into separate areas based on connectivity
    region_face=face(region_face_idx,:);
    region_face_num=size(region_face_idx,1);
    connect_region=zeros(region_face_num,region_face_num+1);
    connect_region_count=0;

    region_face_sign=zeros(region_face_num,1);
    
    while sum(region_face_sign)~=region_face_num
        % Select a seed face for the new connected region and perform the corresponding initialization
        this_face_idx=region_face_idx(find(~region_face_sign,1),1);
        connect_region_count=connect_region_count+1;
        connect_region(connect_region_count,1)=1;
        connect_region(connect_region_count,2)=this_face_idx;
        region_face_sign(find(~region_face_sign,1),1)=1;

        while size(this_face_idx,2)~=0
            new_face_idx=zeros(1,50*size(this_face_idx,1));
            new_face_count=0;
            for i=1:size(this_face_idx,2)
                edge=zeros(3,2);
                edge(1,:)=[face(this_face_idx(1,i),1),face(this_face_idx(1,i),2)];
                edge(2,:)=[face(this_face_idx(1,i),2),face(this_face_idx(1,i),3)];
                edge(3,:)=[face(this_face_idx(1,i),1),face(this_face_idx(1,i),3)];
                for j=1:3
                    third_face_id=findCertainNeighFace(region_face,region_face_idx,edge(j,:),this_face_idx);
                    if third_face_id==0
                        continue;
                    end
                    if region_face_sign(third_face_id,1)==0
                        new_face_count=new_face_count+1;
                        new_face_idx(1,new_face_count)=region_face_idx(third_face_id,1);
                        region_face_sign(third_face_id,1)=1;
                    end
                end
            end
            this_face_idx=new_face_idx(1,1:new_face_count);
            connect_region(connect_region_count,connect_region(connect_region_count,1)+2:...
                connect_region(connect_region_count,1)+1+new_face_count)=this_face_idx;
            connect_region(connect_region_count,1)=connect_region(connect_region_count,1)+new_face_count;
        end
    end

    connect_region=connect_region(1:connect_region_count,:);
end