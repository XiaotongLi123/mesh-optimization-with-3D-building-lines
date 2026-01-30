function [opt_face,opt_val]=DP_reconstruct(point,node,neigh_face_normal,is_constrained,weight)
    n=size(point,1);
    opt_face=[];
    opt_val=inf;

    [is_exist_point_1,exist_point_id_1]=isExistPoint(point,node(1,:));
    [is_exist_point_2,exist_point_id_2]=isExistPoint(point,node(end,:));
    line=[node(1,:);node(end,:)];

    if is_exist_point_1==1&&is_exist_point_2==1
        left=exist_point_id_1(1,1);
        right=exist_point_id_2(1,1);
        node=node(2:end-1,:);
        [opt_face,opt_val]=update(point,node,neigh_face_normal,is_constrained,line,left,right,opt_face,opt_val,weight,1,1);
        return;
    end

    KDmodel=KDTreeSearcher(point);
    neigh=min(3,n);

    if is_exist_point_1==1&&is_exist_point_2==0
        left=exist_point_id_1(1,1);
        index=knnsearch(KDmodel,node(end,:),'K',neigh);
        node=node(2:end,:);
        for i=1:neigh
            right=index(1,i);
            if right==left
                continue;
            end
            [opt_face,opt_val]=update(point,node,neigh_face_normal,is_constrained,line,left,right,opt_face,opt_val,weight,1,0);
            if opt_val~=Inf
                break;
            end
        end
        return;
    end

    if is_exist_point_2==1&&is_exist_point_1==0
        right=exist_point_id_2(1,1);
        index=knnsearch(KDmodel,node(1,:),'K',neigh);
        node=node(1:end-1,:);

        for i=1:neigh
            left=index(1,i);
            if left==right
                continue;
            end
            [opt_face,opt_val]=update(point,node,neigh_face_normal,is_constrained,line,left,right,opt_face,opt_val,weight,0,1);
            if opt_val~=Inf
                break;
            end
        end
        return;
    end 

    index_1=knnsearch(KDmodel,node(1,:),'K',neigh);
    index_2=knnsearch(KDmodel,node(end,:),'K',neigh);
    for i=1:neigh
        for j=1:neigh
            left=index_1(1,i);
            right=index_2(1,j);
            if left==right
                continue;
            end
            [opt_face,opt_val]=update(point,node,neigh_face_normal,is_constrained,line,left,right,opt_face,opt_val,weight,0,0);
            if opt_val~=Inf
                break;
            end
        end
    end
end

function [opt_face,opt_val]=update(point,node,neigh_face_normal,is_constrained,line,left,right,opt_face,opt_val,weight,left_sign,right_sign)
    pts_num=size(point,1);
    node_num=size(node,1);
    n=pts_num;

    if left<right
        point_1=[point(left:right,:);lineDirectionChange(node)];
        neigh_face_normal_1=zeros(size(point_1,1),3);
        neigh_face_normal_1(1:right-left,:)=neigh_face_normal(left:right-1,:);
        is_constrained_1=zeros(size(point_1,1),1);
        is_constrained_1(1:right-left,:)=is_constrained(left:right-1,:);
        point_2=[point(right:pts_num,:);point(1:left,:);node];
        neigh_face_normal_2=zeros(size(point_2,1),3);
        neigh_face_normal_2(1:pts_num-right+left,:)=[neigh_face_normal(right:pts_num,:);neigh_face_normal(1:left-1,:)];
        is_constrained_2=zeros(size(point_2,1),1);
        is_constrained_2(1:pts_num-right+left,:)=[is_constrained(right:pts_num,:);is_constrained(1:left-1,:)];
        map_1=[(left:right)';(n+node_num:-1:n+1)'];
        map_2=[(right:pts_num)';(1:left)';(n+1:n+node_num)'];
        [~,distribution_1]=ransacGetPlaneNormal(point_1,[line(end,:);line(1,:)],1);
        [~,distribution_2]=ransacGetPlaneNormal(point_2,[line(1,:);line(end,:)],2);
    else
        point_1=[point(right:left,:);node];
        neigh_face_normal_1=zeros(size(point_1,1),3);
        neigh_face_normal_1(1:left-right,:)=neigh_face_normal(right:left-1,:);
        is_constrained_1=zeros(size(point_1,1),1);
        is_constrained_1(1:left-right,:)=is_constrained(right:left-1,:);
        point_2=[point(left:pts_num,:);point(1:right,:);lineDirectionChange(node)];
        neigh_face_normal_2=zeros(size(point_2,1),3);
        neigh_face_normal_2(1:pts_num-left+right,:)=[neigh_face_normal(left:pts_num,:);neigh_face_normal(1:right-1,:)];
        is_constrained_2=zeros(size(point_2,1),1);
        is_constrained_2(1:pts_num-left+right,:)=[is_constrained(left:pts_num,:);is_constrained(1:right-1,:)];
        map_1=[(right:left)';(n+1:n+node_num)'];
        map_2=[(left:pts_num)';(1:right)';(n+node_num:-1:n+1)'];
        [~,distribution_1]=ransacGetPlaneNormal(point_1,[line(1,:);line(end,:)],1);
        [~,distribution_2]=ransacGetPlaneNormal(point_2,[line(end,:);line(1,:)],2);
    end
 
    new_neigh_face_normal_1=neigh_face_normal_1;
    new_is_constrained_1=is_constrained_1;
    new_neigh_face_normal_2=neigh_face_normal_2;
    new_is_constrained_2=is_constrained_2;
    
    [opt_face_1,opt_val_1]=DP_ploygonReconstruct(point_1,neigh_face_normal_1,is_constrained_1,1,weight,distribution_1);
    opt_face_1_c=facepointUpdate(opt_face_1,map_1);
    
    if left<right
        for i=pts_num-right+left+1:size(point_2,1)
            if i~=size(point_2,1)
                edge_neigh_face_id=findEdgeNeighFace(opt_face_1_c,[map_2(i,1),map_2(i+1,1)]);
                if i~=pts_num-right+left+1
                    new_is_constrained_2(i,1)=1;
                else
                    if left_sign==1
                        new_is_constrained_2(i,1)=1;
                    end
                end
            else
                edge_neigh_face_id=findEdgeNeighFace(opt_face_1_c,[map_2(i,1),map_2(1,1)]);
                if right_sign==1
                    new_is_constrained_2(i,1)=1;
                end
            end    
            if size(edge_neigh_face_id,1)==0
                continue;
            end
            new_neigh_face_normal_2(i,:)=getFacesNormal([point;node],opt_face_1_c(edge_neigh_face_id,:));
        end
    else
        for i=pts_num-left+right+1:size(point_2,1)
            if i~=size(point_2,1)
                edge_neigh_face_id=findEdgeNeighFace(opt_face_1_c,[map_2(i,1),map_2(i+1,1)]);
                if i~=pts_num-left+right+1
                    new_is_constrained_2(i,1)=1;
                else
                    if right_sign==1
                        new_is_constrained_2(i,1)=1;
                    end
                end
            else
                edge_neigh_face_id=findEdgeNeighFace(opt_face_1_c,[map_2(i,1),map_2(1,1)]);
                if left_sign==1
                    new_is_constrained_2(i,1)=1;
                end
            end    
            if size(edge_neigh_face_id,1)==0
                continue;
            end
            new_neigh_face_normal_2(i,:)=getFacesNormal([point;node],opt_face_1_c(edge_neigh_face_id,:));
        end
    end
    [opt_face_2,opt_val_2]=DP_ploygonReconstruct(point_2,new_neigh_face_normal_2,new_is_constrained_2,1,weight,distribution_2);
    opt_face_2_c=facepointUpdate(opt_face_2,map_2);

    this_opt_val=opt_val_1+opt_val_2;

    opt_face_c=[opt_face_1_c;opt_face_2_c];
    normal_angle=regionFaceNormalAngle([point;node],opt_face_c,(1:size(opt_face_c,1))');

    if max(normal_angle)>75
        [opt_face_2,opt_val_2]=DP_ploygonReconstruct(point_2,neigh_face_normal_2,is_constrained_2,1,weight,distribution_2);
        opt_face_2_c=facepointUpdate(opt_face_2,map_2);
        
        if left<right
            for i=right-left+1:size(point_1,1)
                if i~=size(point_1,1)
                    edge_neigh_face_id=findEdgeNeighFace(opt_face_2_c,[map_1(i,1),map_1(i+1,1)]);
                    if i~=right-left+1
                        new_is_constrained_1(i,1)=1;
                    else
                        if right_sign==1
                            new_is_constrained_1(i,1)=1;
                        end
                    end
                else
                    edge_neigh_face_id=findEdgeNeighFace(opt_face_2_c,[map_1(i,1),map_1(1,1)]);
                    if left_sign==1
                        new_is_constrained_1(i,1)=1;
                    end
                end    
                if size(edge_neigh_face_id,1)==0
                    continue;
                end
                new_neigh_face_normal_1(i,:)=getFacesNormal([point;node],opt_face_2_c(edge_neigh_face_id,:));
            end
        else
            for i=left-right+1:size(point_1,1)
                if i~=size(point_1,1)
                    edge_neigh_face_id=findEdgeNeighFace(opt_face_2_c,[map_1(i,1),map_1(i+1,1)]);
                    if i~=left-right+1
                        new_is_constrained_1(i,1)=1;
                    else
                        if left_sign==1
                            new_is_constrained_1(i,1)=1;
                        end
                    end
                else
                    edge_neigh_face_id=findEdgeNeighFace(opt_face_2_c,[map_1(i,1),map_1(1,1)]);
                    if right_sign==1
                        new_is_constrained_1(i,1)=1;
                    end
                end    
                if size(edge_neigh_face_id,1)==0
                    continue;
                end
                new_neigh_face_normal_1(i,:)=getFacesNormal([point;node],opt_face_2_c(edge_neigh_face_id,:));
            end
        end
        [opt_face_1,opt_val_1]=DP_ploygonReconstruct(point_1,new_neigh_face_normal_1,new_is_constrained_1,1,weight,distribution_1);
        opt_face_1_c=facepointUpdate(opt_face_1,map_1);
    
        new_this_opt_val=opt_val_1+opt_val_2;
        if new_this_opt_val<this_opt_val
            opt_face_c=[opt_face_1_c;opt_face_2_c];
            this_opt_val=new_this_opt_val;
        end
    end

    if this_opt_val<opt_val
        opt_val=this_opt_val;
        opt_face=opt_face_c;
    end
end

function [node_opposite]=lineDirectionChange(node)
    node_num=size(node,1);
    node_opposite=zeros(node_num,3);
    for i=1:node_num
        node_opposite(i,:)=node(node_num+1-i,:);
    end
end

% Initialize the distribution of dihedral angles and assess the flatness of the LSM
function [normal,angle_d]=ransacGetPlaneNormal(point,line)
    pts_num=size(point,1);
    vector=zeros(pts_num,3);

    for i=1:pts_num
        [vector(i,:)]=point2lineVec(point(i,:),line);
    end
    
    center=ransacFindCenter(vector);
    
    normal=crossProduct([line(2,:)-line(1,:);center]);
    if norm(normal)~=0
        normal=normal/norm(normal);
    end
    angle_d=zeros(pts_num,1);
    count=0;
    for i=1:pts_num
        if norm(vector(i,:))~=0
            if getVecAngle(center,vector(i,:))*180/pi<60
                count=count+1;
                angle_d(count,1)=getVecAngle(center,vector(i,:))*180/pi;
            end
        end
    end
    angle_d=angle_d(1:count,1);
end

% Consider both planar and undulating cases
function [center]=ransacFindCenter(data)
    sigma=7.5;
    data_num=size(data,1);
    sign=zeros(data_num,1);
    for i=1:data_num
        if any(data(i,:)~=0)
            sign(i,1)=1;
        end
    end
    data=data(sign==1,:);

    center=mean(data,1);
    data_num=size(data,1);
    pretotal=round(data_num*0.6)-0.1;
    count=0;
    for j=1:data_num
        angle=zeros(data_num,1);
        for m=1:data_num
            angle(m,1)=getVecAngle(data(j,:),data(m,:))*180/pi;
        end
        total=sum(angle<=sigma);
        if total>pretotal
            pretotal=total;
            center=data(j,:);
            count=count+1;
        end
        if total==pretotal
            for m=1:3
                center(1,m)=getMean(center(1,m),count,data(j,m));
                count=count+1;
            end
        end
    end
end

function [vector]=point2lineVec(point,line)
    if all(point==line(1,:))||all(point==line(2,:))
        vector=zeros(1,3);
        return;
    end
    [a,b,c]=LineEquation(line(1,:),line(2,:));

    d=-(a*point(1,1)+b*point(1,2)+c*point(1,3));

    t=-(a*line(1,1)+b*line(1,2)+c*line(1,3)+d)/(a^2+b^2+c^2);
    O=[line(1,1)+a*t,line(1,2)+b*t,line(1,3)+c*t];

    vector=point-O;
    scale=sqrt(vector(1,1)^2+vector(1,2)^2+vector(1,3)^2);
    if scale~=0
        for j=1:3
            vector(1,j)=vector(1,j)/scale; 
        end
    end
end

function [new_mean]=getMean(mean,n,new_data)
    new_mean=(n*mean+new_data)/(n+1);
end