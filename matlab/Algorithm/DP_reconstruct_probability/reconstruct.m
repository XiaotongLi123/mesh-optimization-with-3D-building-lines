function [new_face,new_point,normal_angle,constrained_line,error]=reconstruct(input_face,input_point,input_region_face_idx,line,PFneighbor,constrained_line,interval)   
    input_face_num=size(input_face,1);
    input_region_face_num=size(input_region_face_idx,1);

    new_face_ori=zeros(input_face_num+10*input_region_face_num,3);
    new_face_ori(1:input_face_num,:)=input_face;
    error=0;

    region_face=input_face(input_region_face_idx,:);
    [region_point,~]=getRegionPE(region_face);
    [edge_point_sort,edge_sort]=getOutlineSort(region_face);
    outline_num=size(edge_point_sort,1);
    
    [constrained_outline]=findConstrainedOutline(input_point,edge_point_sort,line);
    if size(constrained_outline,1)~=0
        new_face=input_face;
        new_point=input_point;
        constrained_line=[constrained_line;constrained_outline];
        normal_angle=zeros(1,3);
        return;
    end
    

    [weight_1,weight_2]=getWeight(input_point,input_face,input_region_face_idx,0.5,0.5);
    weight=weight_2/weight_1;
    
    enlarge_face_idx=getPointNeighFace(region_point(:,1),PFneighbor);
    neigh_face_normal=zeros(outline_num,3);
    is_constrained=zeros(outline_num,1);
    for i=1:outline_num
        edge_neigh_face_id=findEdgeNeighFace(input_face(enlarge_face_idx,:),edge_sort(i,:));
        edge_neigh_face_idx=enlarge_face_idx(edge_neigh_face_id,:);
        if size(edge_neigh_face_idx,1)==1
            continue;
        end
        neigh_face_idx=setdiff(edge_neigh_face_idx,input_region_face_idx);
        neigh_face=input_face(neigh_face_idx,:);
        is_constrained(i,1)=isExistLine(constrained_line,edge_sort(i,:));
        neigh_face_normal(i,:)=getFacesNormal(input_point,neigh_face);
    end

    node=getNode(line,interval);

    [opt_face,opt_val]=DP_reconstruct(input_point(edge_point_sort,:),node,neigh_face_normal,is_constrained,weight);
    if opt_val==inf
        error=1;
        new_face=input_face;
        new_point=input_point;
        normal_angle=zeros(1,3);
        return;
    end

    [map,new_constrained_line,new_point]=postProcess(input_point,node,edge_point_sort);
    constrained_line=[constrained_line;new_constrained_line];
    opt_face_c=faceChange(opt_face,map);
    reconstruct_face_num=size(opt_face_c,1);

    new_face_ori(input_region_face_idx,:)=0;
    new_face_num=input_face_num-input_region_face_num+reconstruct_face_num;
    new_face=zeros(new_face_num,3);
    count=0;
    for i=1:input_face_num
        if new_face_ori(i,1)~=0
            count=count+1;
            new_face(count,:)=new_face_ori(i,:);
        end
    end
    new_face(count+1:new_face_num,:)=opt_face_c;

    for i=1:size(region_point,1)
        if region_point(i,2)==1
            new_point(region_point(i,1),:)=0;
        end
    end

    normal_angle=regionFaceNormalAngle(new_point,new_face,(new_face_num-reconstruct_face_num+1:new_face_num)');
end

function [new_face]=faceChange(old_face,map)
    face_num=size(old_face,1);
    new_face=zeros(face_num,3);
    for i=1:face_num
        for j=1:3
            new_face(i,j)=map(old_face(i,j),1);
        end
    end
end

function [map,constrained_line,new_point]=postProcess(input_point,node,edge_point_sort)
    [is_exist_point_1,exist_point_id_1]=isExistPoint(input_point(edge_point_sort,:),node(1,:));
    [is_exist_point_2,exist_point_id_2]=isExistPoint(input_point(edge_point_sort,:),node(end,:));
    node_num=size(node,1);
    input_point_num=size(input_point,1);

    constrained_line=zeros(node_num-1,2);

    if is_exist_point_1==0
        if is_exist_point_2==0
            map=[edge_point_sort;(input_point_num+1:input_point_num+node_num)'];
            new_point=[input_point;node];

            for i=1:node_num-1
                constrained_line(i,:)=[input_point_num+i,input_point_num+i+1];
            end
        else
            map=[edge_point_sort;(input_point_num+1:input_point_num+node_num-1)'];
            new_point=[input_point;node(1:end-1,:)];

            for i=1:node_num-2
                constrained_line(i,:)=[input_point_num+i,input_point_num+i+1];
            end
            constrained_line(node_num-1,:)=[input_point_num+node_num-1,edge_point_sort(exist_point_id_2,1)];
        end
    else
        if is_exist_point_2==0
            map=[edge_point_sort;(input_point_num+1:input_point_num+node_num-1)'];
            new_point=[input_point;node(2:end,:)];

            constrained_line(1,:)=[edge_point_sort(exist_point_id_1,1),input_point_num+1];
            for i=1:node_num-2
                constrained_line(i+1,:)=[input_point_num+i,input_point_num+i+1];
            end
        else
            map=[edge_point_sort;(input_point_num+1:input_point_num+node_num-2)'];
            new_point=[input_point;node(2:end-1,:)];
            if node_num>2
                constrained_line(1,:)=[edge_point_sort(exist_point_id_1,1),input_point_num+1];
                for i=1:node_num-3
                    constrained_line(i+1,:)=[input_point_num+i,input_point_num+i+1];
                end
                constrained_line(node_num-1,:)=[input_point_num+node_num-2,edge_point_sort(exist_point_id_2,1)];
            else
                constrained_line(1,:)=[edge_point_sort(exist_point_id_1,1),edge_point_sort(exist_point_id_2,1)];
            end
        end
    end
end

function [constrained_outline]=findConstrainedOutline(point,outline_point_sort,line)
    outline_point_num=size(outline_point_sort,1);
    constrained_outline=[];
    left=0;
    right=0;
    for i=1:outline_point_num
        if point2pointDist(point(outline_point_sort(i,1),:),line(1,:))<1e-4
            left=i;
        end
        if point2pointDist(point(outline_point_sort(i,1),:),line(2,:))<1e-4
            right=i;
        end
    end
    overlapping_sign=ones(2,1);
    if left~=right&&left~=0&&right~=0
        for i=1:outline_point_num
            if i>=min(left,right)&&i<=max(left,right)
                if point2segmentDist(point(outline_point_sort(i,1),:),line)>1e-4
                    overlapping_sign(1,1)=0;
                end
            else
                if point2segmentDist(point(outline_point_sort(i,1),:),line)>1e-4
                    overlapping_sign(2,1)=0;
                end
            end
        end
    end
    if overlapping_sign(1,1)==1&&overlapping_sign(2,1)==0
        constrained_outline=zeros(max(left,right)-min(left,right),2);
        for i=min(left,right):max(left,right)-1
            constrained_outline(i+1-min(left,right),:)=[outline_point_sort(i,1),outline_point_sort(i+1,1)];
        end
    end
    if overlapping_sign(2,1)==1&&overlapping_sign(1,1)==0
        constrained_outline=zeros(outline_point_num-(max(left,right)-min(left,right)),2);
        for i=1:min(left,right)-1
            constrained_outline(i,:)=[outline_point_sort(i,1),outline_point_sort(i+1,1)];
        end
        for i=max(left,right):outline_point_num-1
            constrained_outline(i+min(left,right)-1,:)=[outline_point_sort(i,1),outline_point_sort(i+1,1)];
        end
        constrained_outline(outline_point_num-(max(left,right)-min(left,right)),:)=[outline_point_sort(outline_point_num,1),outline_point_sort(1,1)];
    end
end