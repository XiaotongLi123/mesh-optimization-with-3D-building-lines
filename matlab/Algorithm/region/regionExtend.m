function [extend_region_face_idx,error]=regionExtend(face,knn_face_idx,region_face_idx,normal,dist,PFneighbor,threshold_1,threshold_2)
    % LSM growing
    % dist contains the distance from each face in initialized LSM to the constrained line segments
    sign=1; % Exit the loop when sign = 0
    knn_face=face(knn_face_idx,:);
    
    error=0;
    loop=0;
    region_face_num=size(region_face_idx,1);
    region_face_idx=[region_face_idx;zeros(size(face,1),1)];
    
    while sign==1&&loop<1000
        
        region_face=face(region_face_idx(1:region_face_num,1),:);   
        [region_point,~]=getRegionPE(region_face);
     
        enlarge_face_idx=getPointNeighFace(region_point(:,1),PFneighbor); % Region growth
        enlarge_face_idx=intersect(enlarge_face_idx,knn_face_idx);
    
        growth_face_idx=setdiff(enlarge_face_idx,region_face_idx(1:region_face_num,1)); % Take the set difference between two sets to obtain the IDs of the grown faces
        growth_face=face(growth_face_idx,:);

        % Determine whether each grown face can be further expanded in the region
        growth_face_num=size(growth_face_idx,1);
        isExtend=zeros(growth_face_num,1);
        for i=1:growth_face_num
            % For each grown face, find all neighboring faces within the region through shared edges
            edge=[growth_face(i,2),growth_face(i,3);
                growth_face(i,3),growth_face(i,1);
                growth_face(i,1),growth_face(i,2)];
            neigh_region_face_idx=zeros(100000,1);
            count=0;
            for j=1:3
                [neigh_face_id]=findEdgeNeighFace(knn_face,edge(j,:));
                neigh_face_idx=knn_face_idx(neigh_face_id,1);
                neigh_face_idx=setdiff(neigh_face_idx,growth_face_idx(i,1));% Remove the current face
                try
                    neigh_region_face_idx(count+1:count+size(intersect(neigh_face_idx,region_face_idx(1:region_face_num,1)),1),1)=intersect(neigh_face_idx,region_face_idx(1:region_face_num,1));
                catch
                    neigh_region_face_idx(count+1:count+size(intersect(neigh_face_idx,region_face_idx(1:region_face_num,1)),2),1)=intersect(neigh_face_idx,region_face_idx(1:region_face_num,1));
                end
                count=count+size(intersect(neigh_face_idx,region_face_idx(1:region_face_num,1)),1);
            end
            neigh_region_face_idx=neigh_region_face_idx(1:count,1);

            % Determine whether the current face should grow based on the distance from the constrained lines 
            % to the current face and its neighboring faces in the influence region
            if size(neigh_region_face_idx,1)==0||size(neigh_region_face_idx,2)==0
                continue;
            end
            angle=zeros(size(neigh_region_face_idx,1),1);
            for j=1:size(neigh_region_face_idx,1)
                angle(j,1)=getVecAngle(normal(knn_face_idx==neigh_region_face_idx(j,1),:),normal(knn_face_idx==growth_face_idx(i,1),:))*180/pi;
            end
            min_angle=min(angle);
            if min_angle<threshold_1&&dist(knn_face_idx==growth_face_idx(i,1),1)<threshold_2
                isExtend(i,1)=1;
            end
        end
        extend_face_num=sum(isExtend);

        if extend_face_num==0
            sign=0;
        end
        region_face_idx(region_face_num+1:region_face_num+extend_face_num,1)=growth_face_idx(isExtend==1,1);
        region_face_num=region_face_num+extend_face_num;

        loop=loop+1;
    end

    if loop==1000
        error=1;
    end
    extend_region_face_idx=region_face_idx(1:region_face_num,1);
end