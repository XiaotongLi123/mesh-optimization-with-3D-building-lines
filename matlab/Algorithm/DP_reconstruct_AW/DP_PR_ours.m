function [opt_face,opt_val]=DP_PR_ours(point,neigh_face_normal,is_constrained,weight_1,weight_2)
    % Dynamic programming polygon triangulation based on global optimization of triangle perimeter and dihedral angle sum
    % Inputs: boundary vertices of the polygon, initial normals, and two weight factors
    % Outputs: reconstructed triangle faces and the corresponding objective function value
    pts_num=size(point,1);
    n=pts_num;
    face=cell((n-1)*(n-2)/2,1); % A convex polygon with n edges contains n-2 triangles
    % Pre-allocate memory for matrices
    for s=3:n
        for i=1:n-s+1
            face{getStoragePos(i,s,n),1}=zeros(s-2,s-2,3);
        end
    end
    
    weight=weight_2/weight_1;
    %----------------Initialization----------------
    % Build a list of distances between all pairs of points
    dist=zeros(n,n);
    for i=1:n
        for j=1:n
            dist(i,j)=point2pointDist(point(i,:),point(j,:));
        end
    end
    % Build triangle areas and corresponding normal vectors, ensuring collinear points do not form a triangle
    tri_normal=zeros(n^3,3);
    is_tri=zeros(n^3,1);
    for i=1:n
        for j=1:n
            for k=1:n
                if i==j||j==k||k==i
                    continue;
                end
                angle_1=getLineAngle(point(j,:)-point(i,:),point(i,:)-point(k,:));
                angle_2=getLineAngle(point(k,:)-point(j,:),point(j,:)-point(i,:));
                angle_3=getLineAngle(point(i,:)-point(k,:),point(k,:)-point(j,:));
                if min([angle_1;angle_2;angle_3])>0.1/180*pi                   
                    is_tri((i-1)*n^2+(j-1)*n+k,1)=1;
                end
                if is_tri((i-1)*n^2+(j-1)*n+k,1)==1
                    tri_normal((i-1)*n^2+(j-1)*n+k,:)=getFacesNormal([point(i,:);point(j,:);point(k,:)],[1,2,3]);
                    if norm(tri_normal((i-1)*n^2+(j-1)*n+k,:))~=0
                        tri_normal((i-1)*n^2+(j-1)*n+k,:)=tri_normal((i-1)*n^2+(j-1)*n+k,:)/norm(tri_normal((i-1)*n^2+(j-1)*n+k,:));
                    end
                end    
            end
        end
    end
    %----------------Initialize the DP subproblem corresponding to a 3-edge polygon----------------
    opt_func=zeros(n,n,n-2);
    opt_length=zeros(n,n,n-2);
    opt_Dangle=zeros(n,n,n-2);
    neigh_normal=zeros(n^2,n-2,3);
    for i=1:n-2
        % If three points are collinear, a triangle cannot be constructed
        if is_tri((i-1)*n^2+i*n+i+2,1)==0
            opt_func(i,3,1)=Inf;
            continue;
        end
        face_temp=zeros(1,1,3);
        face_temp(1,1,:)=[i,i+1,i+2];
        face{getStoragePos(i,3,n),1}=face_temp;
        neigh_angle_1=0;
        neigh_angle_2=0;
        if ~all(neigh_face_normal(i,:)==0)
            neigh_angle_1=getVecAngle(tri_normal((i-1)*n^2+i*n+i+2,:),neigh_face_normal(i,:))*180/pi;
            neigh_angle_1=punishFun(neigh_angle_1,is_constrained(i,1),weight);
        end
        if ~all(neigh_face_normal(i+1,:)==0)
            neigh_angle_2=getVecAngle(tri_normal((i-1)*n^2+i*n+i+2,:),neigh_face_normal(i+1,:))*180/pi;
            neigh_angle_2=punishFun(neigh_angle_2,is_constrained(i+1,1),weight);
        end
        opt_length(i,3,1)=dist(i,i+1)+dist(i+1,i+2)+dist(i+2,i);
        opt_func(i,3,1)=weight_1*(dist(i,i+1)+dist(i+1,i+2)+dist(i+2,i))+weight_2*(neigh_angle_1+neigh_angle_2);
        neigh_normal((i-1)*n+3,1,:)=tri_normal((i-1)*n^2+i*n+i+2,:);
    end

    %----------------Set up the dynamic programming recurrence for solving----------------
    for s=4:n
        perimeter=zeros(10000,1);
        Dangle=zeros(10000,1);
        this_count=0;
        for i=1:n-s+1
            face_temp=zeros(s-2,s-2,3);
            for k=1:s-2
                % Exclude the case where the newly formed triangle is collinear
                if is_tri((i-1)*n^2+(i+k-1)*n+i+s-1,1)==0
                    opt_func(i,s,k)=Inf;
                    continue;
                end
                if k~=1&&k~=s-2
                    temp=zeros(((k+1)-2)*(s-(k+1)-1),1);
                    temp_length=zeros(((k+1)-2)*(s-(k+1)-1),1);
                    temp_Dangle=zeros(((k+1)-2)*(s-(k+1)-1),1);
                    for p=1:(k+1)-2
                        for q=1:s-(k+1)-1
                            left=zeros(1,3);
                            right=zeros(1,3);
                            for a=1:3
                                left(1,a)=neigh_normal((i-1)*n+k+1,p,a);
                                right(1,a)=neigh_normal((i+k-1)*n+s-k,q,a);
                            end
                            if any(left~=0)&&any(right~=0)
                                temp_length((p-1)*(s-(k+1)-1)+q,1)=opt_length(i,k+1,p)+opt_length(i+k,s-k,q)+(dist(i,i+k)+dist(i+k,i+s-1)+dist(i+s-1,i));
                                temp_Dangle((p-1)*(s-(k+1)-1)+q,1)=opt_Dangle(i,k+1,p)+opt_Dangle(i+k,s-k,q)+getVecAngle(tri_normal((i-1)*n^2+(i+k-1)*n+i+s-1,:),left)*180/pi+...
                                    +getVecAngle(tri_normal((i-1)*n^2+(i+k-1)*n+i+s-1,:),right)*180/pi;
                                temp((p-1)*(s-(k+1)-1)+q,1)=opt_func(i,k+1,p)+opt_func(i+k,s-k,q)+(dist(i,i+k)+...
                                    dist(i+k,i+s-1)+dist(i+s-1,i))+weight*punishFun(getVecAngle(tri_normal((i-1)*n^2+(i+k-1)*n+i+s-1,:),left)*180/pi,-1,weight)+...
                                    +weight*punishFun(getVecAngle(tri_normal((i-1)*n^2+(i+k-1)*n+i+s-1,:),right)*180/pi,-1,weight);
                            else
                                temp((p-1)*(s-(k+1)-1)+q,1)=Inf;
                            end
                        end
                    end
                    [opt_func(i,s,k),min_id]=min(temp);
                    opt_length(i,s,k)=temp_length(min_id,1);
                    opt_Dangle(i,s,k)=temp_Dangle(min_id,1);
                    if opt_func(i,s,k)~=Inf
                        this_count=this_count+1;
                        perimeter(this_count,1)=opt_length(i,s,k);
                        Dangle(this_count,1)=opt_Dangle(i,s,k);
                    end
        
                    [r,c]=getRowCol(s-(k+1)-1,min_id);
                    face_left=face{getStoragePos(i,k+1,n),1};
                    face_right=face{getStoragePos(i+k,s-k,n),1};
                    face_temp(k,1:(k+1-2),:)=face_left(r,:,:); % Triangle faces corresponding to the first subproblem
                    face_temp(k,(k+1-2+1):s-2-1,:)=face_right(c,:,:); % Triangle faces corresponding to the second subproblem
                    face_temp(k,s-2,:)=[i,i+k,i+s-1]; % Construct a triangle using points i, i+k, and i+s-1
                else
                    if k==1
                        temp=zeros(s-(k+1)-1,1);
                        temp_length=zeros(s-(k+1)-1,1);
                        temp_Dangle=zeros(s-(k+1)-1,1);
                        for q=1:s-(k+1)-1
                            right=zeros(1,3);
                            for a=1:3
                                right(1,a)=neigh_normal((i+k-1)*n+s-k,q,a);
                            end
                            if any(right~=0)
                                temp_length(q,1)=opt_length(i+k,s-k,q)+(dist(i,i+k)+dist(i+k,i+s-1)+dist(i+s-1,i));
                                temp_Dangle(q,1)=opt_Dangle(i+k,s-k,q)+getVecAngle(tri_normal((i-1)*n^2+(i+k-1)*n+i+s-1,:),right)*180/pi;
                                temp(q,1)=opt_func(i+k,s-k,q)+(dist(i,i+k)+dist(i+k,i+s-1)+dist(i+s-1,i))+...
                                    weight*punishFun(getVecAngle(tri_normal((i-1)*n^2+(i+k-1)*n+i+s-1,:),right)*180/pi,-1,weight);
                            else
                                temp(q,1)=Inf;
                                continue;
                            end
                            if any(neigh_face_normal(i,:)~=0)
                                temp(q,1)=temp(q,1)+weight*punishFun(getVecAngle(tri_normal((i-1)*n^2+(i+k-1)*n+i+s-1,:),neigh_face_normal(i,:))*180/pi,is_constrained(i,1),weight);
                            end
                        end
                        [opt_func(i,s,k),min_id]=min(temp);
                        opt_length(i,s,k)=temp_length(min_id,1);
                        opt_Dangle(i,s,k)=temp_Dangle(min_id,1);
                        if opt_func(i,s,k)~=Inf
                            this_count=this_count+1;
                            perimeter(this_count,1)=opt_length(i,s,k);
                            Dangle(this_count,1)=opt_Dangle(i,s,k);
                        end

                        face_right=face{getStoragePos(i+k,s-k,n),1};
                        face_temp(k,(k+1-2+1):s-2-1,:)=face_right(min_id,:,:); % Triangle faces corresponding to the second subproblem
                        face_temp(k,s-2,:)=[i,i+k,i+s-1]; % Construct a triangle using points i, i+k, and i+s-1
                    else
                        temp=zeros((k+1)-2,1);
                        temp_length=zeros((k+1)-2,1);
                        temp_Dangle=zeros((k+1)-2,1);
                        for p=1:(k+1)-2
                            left=zeros(1,3);
                            for a=1:3
                                left(1,a)=neigh_normal((i-1)*n+k+1,p,a);   
                            end
                            if any(left~=0)
                                temp_length(p,1)=opt_length(i,k+1,p)+(dist(i,i+k)+dist(i+k,i+s-1)+dist(i+s-1,i));
                                temp_Dangle(p,1)=opt_Dangle(i,k+1,p)+getVecAngle(tri_normal((i-1)*n^2+(i+k-1)*n+i+s-1,:),left)*180/pi;
                                temp(p,1)=opt_func(i,k+1,p)+(dist(i,i+k)+dist(i+k,i+s-1)+dist(i+s-1,i))+...
                                    weight*punishFun(getVecAngle(tri_normal((i-1)*n^2+(i+k-1)*n+i+s-1,:),left)*180/pi,-1,weight);
                            else
                                temp(p,1)=Inf;
                                continue;
                            end
                            if any(neigh_face_normal(i+s-2,:)~=0) 
                                temp(p,1)=temp(p,1)+weight*punishFun(getVecAngle(tri_normal((i-1)*n^2+(i+k-1)*n+i+s-1,:),neigh_face_normal(i+s-2,:))*180/pi,is_constrained(i+s-2,1),weight);
                            end
                        end
                        [opt_func(i,s,k),min_id]=min(temp);
                        opt_length(i,s,k)=temp_length(min_id,1);
                        opt_Dangle(i,s,k)=temp_Dangle(min_id,1);
                        if opt_func(i,s,k)~=Inf
                            this_count=this_count+1;
                            perimeter(this_count,1)=opt_length(i,s,k);
                            Dangle(this_count,1)=opt_Dangle(i,s,k);
                        end
                        face_left=face{getStoragePos(i,k+1,n),1};
                        face_temp(k,1:(k+1-2),:)=face_left(min_id,:,:); % Triangle faces corresponding to the first subproblem
                        face_temp(k,s-2,:)=[i,i+k,i+s-1]; % Construct a triangle using points i, i+k, and i+s-1
                    end
                end
                face{getStoragePos(i,s,n),1}=face_temp;
                neigh_normal((i-1)*n+s,k,:)=tri_normal((i-1)*n^2+(i+k-1)*n+i+s-1,:);
                % When s == n, an additional dihedral angle needs to be considered
                if s==n
                    if any(neigh_face_normal(n,:)~=0)
                        opt_func(i,s,k)=opt_func(i,s,k)+weight*punishFun(getVecAngle(tri_normal((i-1)*n^2+(i+k-1)*n+i+s-1,:),neigh_face_normal(n,:))*180/pi,is_constrained(n,1),weight);
                    end
                end
            end            
        end
        weight=1/(sum(Dangle(1:this_count,1)+1)*(s-2)/(sum(perimeter(1:this_count,1))*(s-3)));
    end
    % Return the optimal triangulation score/value
    [opt_val,opt_id]=min(opt_func(1,n,:));
    % Return the face table of the optimal triangulation
    opt_face=zeros(n-2,3);
    opt_face_cell=face{getStoragePos(1,n,n),1};
    for i=1:n-2
        opt_face(i,:)=opt_face_cell(opt_id,i,:);
    end
end

function [angle]=punishFun(angle,sign,weight)
    lambda=1000*pi/weight/180;
    if sign==1
        if angle>135
            angle=lambda*angle;
        else
            angle=0;
        end
    end
    if sign==0
        if angle>105
            angle=lambda*angle;
        else
            angle=0;
        end
    end
    if sign==-1
        if angle>135
            angle=lambda*angle;
        else
            if angle>105
                angle=3*angle;
            end
        end
    end
end

function [row_pos,col_pos]=getRowCol(cols,this_pos)
    col_pos=mod(this_pos,cols);
    if col_pos==0
        col_pos=cols;
    end
    row_pos=(this_pos-col_pos)/cols+1;
end

function [pos]=getStoragePos(i,s,n)
    pos=(s-3)*((n-2)+(n-(s-2)))/2+i;
end

