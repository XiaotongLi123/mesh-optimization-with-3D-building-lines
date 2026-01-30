function [repaired_region_face_idx]=shapeRepair(point,region_face_idx,face,PFneighbor)
    % Shape repair algorithm
    % Perform a single region growing step: check whether all three vertices of a grown triangle are within the current region;
    % if so, include the triangle in the region
    % Repeat the process until no new triangles satisfy the condition
    % Then, iterate over all boundary points and check if the sum of internal angles of adjacent triangles is less than 240 degrees;
    % if not, expand the corresponding vertex
        
    region_face=face(region_face_idx,:);
    [region_point,~]=getRegionPE(region_face);

    % Region growth  
    sign=1;
    loop=0;
    count=size(region_face_idx,1);
    region_face_idx=[region_face_idx;zeros(size(face,1),1)];
    while sign==1&&loop<100
        enlarge_face_idx=getPointNeighFace(region_point(:,1),PFneighbor);

        % Extract the newly grown faces
        grow_face_idx=zeros(size(enlarge_face_idx,1),1);
        grow_face_num=0;
        for i=1:size(enlarge_face_idx,1)
            if ~any(ismember(region_face_idx(1:count,1),enlarge_face_idx(i,1)))
                grow_face_num=grow_face_num+1;
                grow_face_idx(grow_face_num,1)=enlarge_face_idx(i,1);
            end
        end
        grow_face_idx=grow_face_idx(1:grow_face_num,1);

        % Determine whether each grown face satisfies the specified conditions
        isExtend=zeros(grow_face_num,1);
        for i=1:grow_face_num
            this_face=face(grow_face_idx(i,1),:);
            if any(ismember(region_point(:,1),this_face(1,1)))&&...
                    any(ismember(region_point(:,1),this_face(1,2)))&&any(ismember(region_point(:,1),this_face(1,3)))
                isExtend(i,1)=1;
            end
        end

        % Growing
        loop=loop+1;
        if sum(isExtend)==0
            sign=0;
        else
            region_face_idx(count+1:count+sum(isExtend),1)=grow_face_idx(isExtend==1,1);
            count=count+sum(isExtend);
            region_face=face(region_face_idx(1:count),:);
            [region_point,~]=getRegionPE(region_face);
        end
    end
    region_face_idx=region_face_idx(1:count,1);
    % Check the angles at each vertex
    sign=1;
    loop=0;
    while sign==1&&loop<100
        old_region_face_idx=region_face_idx;
        % Identify and extract the boundary points of the given region
        region_face=face(region_face_idx,:);
        [region_point,~]=getRegionPE(region_face);
        outline_point=region_point(region_point(:,2)==0,1);
        outline_point_num=size(outline_point,1);
        for i=1:outline_point_num
            % Extract adjacent faces
            this_neigh_face_idx=PFneighbor(outline_point(i,1),2:PFneighbor(outline_point(i,1),1)+1)';
            % Retrieve all faces inside the region that are adjacent to a given face
            in_this_neigh_face_idx=intersect(this_neigh_face_idx,region_face_idx);
            % Compute the sum of angles
            in_this_neigh_face=face(in_this_neigh_face_idx,:);
            sum_angle=getVertexAngleSum(point,in_this_neigh_face,outline_point(i,1));
            if sum_angle>240
                region_face_idx=union(region_face_idx,this_neigh_face_idx);
            end
        end
        if size(old_region_face_idx,1)==size(region_face_idx,1)
            sign=0;
        end
        loop=loop+1;
    end
    repaired_region_face_idx=region_face_idx;
end

function [sum_angle]=getVertexAngleSum(point,face,this_vertex)
    face_num=size(face,1);
    sum_angle=0;
    for i=1:face_num
        not_vertex=face(i,face(i,:)~=this_vertex)';
        vec_1=point(not_vertex(1,:),:)-point(this_vertex,:);
        vec_2=point(not_vertex(2,:),:)-point(this_vertex,:);
        sum_angle=sum_angle+getVecAngle(vec_1,vec_2);
    end
    sum_angle=sum_angle*180/pi;
end