function [opt_face,opt_val]=DP_ploygonReconstruct(point,neigh_face_normal,is_constrained,weight_1,weight_2,distribution)
    % Perform polygon triangulation using dynamic programming, enhanced with a probabilistic model to reduce search space and improve efficiency
    pts_num=size(point,1);
    n=pts_num;
    face=cell((n-1)*(n-2)/2,1);
    for s=3:n
        for i=1:n-s+1
            face{getStoragePos(i,s,n),1}=zeros(s-2,s-2,3);
        end
    end
    %----------------Initialization----------------
    dist=zeros(n,n);
    for i=1:n
        for j=1:n
            dist(i,j)=point2pointDist(point(i,:),point(j,:));
        end
    end
    tri_normal=zeros(n^3,3);
    is_tri=zeros(n^3,1);

    count=size(distribution,1);
    distribution=[distribution;zeros(999999,1)];
    
    for i=1:n
        for j=i+1:n
            for k=j+1:n
                if i==j||j==k||k==i
                    continue;
                end
                angle_1=getLineAngle(point(j,:)-point(i,:),point(i,:)-point(k,:));
                angle_2=getLineAngle(point(k,:)-point(j,:),point(j,:)-point(i,:));
                angle_3=getLineAngle(point(i,:)-point(k,:),point(k,:)-point(j,:));
                if min([angle_1;angle_2;angle_3])>0.01/180*pi                   
                    is_tri((i-1)*n^2+(j-1)*n+k,1)=1;
                end
                if is_tri((i-1)*n^2+(j-1)*n+k,1)==1
                    tri_normal((i-1)*n^2+(j-1)*n+k,:)=getFacesNormal([point(i,:);point(j,:);point(k,:)],[1,2,3]);
                end    
            end
        end
    end

    %----------------Initialize the DP subproblem corresponding to a 3-edge polygon----------------
    total = 0;          % Number of sub-polygon merges
    back = 0;           % Number of backtracking occurrences
    noback_correct = 0; % Number of times backtracking was rejected and the solution was correct
    noback_10 = 0;      % Number of times backtracking was rejected but the optimal value has ~10% error
    noback_20 = 0;      % Number of times backtracking was rejected but the optimal value has ~20% error
    noback_more = 0;    % Number of times backtracking was rejected but the optimal value error exceeds 20%
    back_noneed = 0;    % Number of unnecessary backtracking occurrences

    opt_func=zeros(n,n,n-2);
    neigh_normal=zeros(n^2,n-2,3);
    neigh_angle=zeros(n^2,n-2,2);
    opt_idx=zeros(n,n);
    for i=1:n-2
        opt_idx(i,3)=1;
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
            neigh_angle_1=punishFun(neigh_angle_1,is_constrained(i,1),weight_2);
        end
        if ~all(neigh_face_normal(i+1,:)==0)
            neigh_angle_2=getVecAngle(tri_normal((i-1)*n^2+i*n+i+2,:),neigh_face_normal(i+1,:))*180/pi;
            neigh_angle_2=punishFun(neigh_angle_2,is_constrained(i+1,1),weight_2);
        end
        opt_func(i,3,1)=weight_1*(dist(i,i+1)+dist(i+1,i+2)+dist(i+2,i))+weight_2*(neigh_angle_1+neigh_angle_2);
        neigh_normal((i-1)*n+3,1,:)=tri_normal((i-1)*n^2+i*n+i+2,:);
    end
    
    %----------------Set up the dynamic programming recurrence for solving----------------
    thre_count=0;
    thre_set=zeros(100,1);
    for s=4:n
        s_count=count;
        THRE=getThreshold(distribution(1:s_count,1),0.05);
        thre_count=thre_count+1;
        thre_set(thre_count,1)=THRE;
        for i=1:n-s+1
            face_temp=zeros(s-2,s-2,3);
            for k=1:s-2
                if is_tri((i-1)*n^2+(i+k-1)*n+i+s-1,1)==0
                    opt_func(i,s,k)=Inf;
                    continue;
                end
                if opt_idx(i,k+1)~=0
                    if opt_func(i,k+1,opt_idx(i,k+1))==Inf
                        opt_func(i,s,k)=Inf;
                        continue;
                    end
                end
                if opt_idx(i+k,s-k)~=0
                    if opt_func(i+k,s-k,opt_idx(i+k,s-k))==Inf
                        opt_func(i,s,k)=Inf;
                        continue;
                    end
                end
                if k~=1&&k~=s-2
                    total=total+2;

                    temp_1=zeros((k+1)-2,1);
                    temp_2=zeros(s-(k+1)-1,1);
                    left=zeros(1,3);
                    right=zeros(1,3);

                    for p=1:(k+1)-2
                        for a=1:3
                            left(1,a)=neigh_normal((i-1)*n+k+1,p,a);
                        end
                        if any(left~=0)
                            temp_1(p,1)=opt_func(i,k+1,p)+weight_2*punishFun(getVecAngle(tri_normal((i-1)*n^2+(i+k-1)*n+i+s-1,:),left)*180/pi,-1,weight_2);
                        else
                            temp_1(p,1)=Inf;
                        end
                    end
                    [opt_val_1,r_1]=min(temp_1);
                    for q=1:s-(k+1)-1
                        for a=1:3
                            right(1,a)=neigh_normal((i+k-1)*n+s-k,q,a);
                        end
                        if any(right~=0) 
                            temp_2(q,1)=opt_func(i+k,s-k,q)+weight_2*punishFun(getVecAngle(tri_normal((i-1)*n^2+(i+k-1)*n+i+s-1,:),right)*180/pi,-1,weight_2);
                        else
                            temp_2(q,1)=Inf;
                        end
                    end
                    [opt_val_2,c_1]=min(temp_2);

                    left_opt_normal=zeros(1,3);
                    right_opt_normal=zeros(1,3);
                    for a=1:3
                        left_opt_normal(1,a)=neigh_normal((i-1)*n+k+1,opt_idx(i,k+1),a);
                        right_opt_normal(1,a)=neigh_normal((i+k-1)*n+s-k,opt_idx(i+k,s-k),a);
                    end

                    opt_ref_1=opt_func(i,k+1,opt_idx(i,k+1))+weight_2*punishFun(getVecAngle(tri_normal((i-1)*n^2+(i+k-1)*n+i+s-1,:),left_opt_normal)*180/pi,-1,weight_2);
                    opt_ref_2=opt_func(i+k,s-k,opt_idx(i+k,s-k))+weight_2*punishFun(getVecAngle(tri_normal((i-1)*n^2+(i+k-1)*n+i+s-1,:),right_opt_normal)*180/pi,-1,weight_2);

                    angle_left=getVecAngle(tri_normal((i-1)*n^2+(i+k-1)*n+i+s-1,:),left_opt_normal)*180/pi;
                    angle_right=getVecAngle(tri_normal((i-1)*n^2+(i+k-1)*n+i+s-1,:),right_opt_normal)*180/pi;
                    if angle_left>=THRE&&angle_right>=THRE
                        opt_func(i,s,k)=opt_val_1+opt_val_2+weight_1*(dist(i,i+k)+dist(i+k,i+s-1)+dist(i+s-1,i));

                        back=back+2;

                        if opt_idx(i,k+1)==r_1
                            back_noneed=back_noneed+1;
                        end
                        if opt_idx(i+k,s-k)==c_1
                            back_noneed=back_noneed+1;
                        end
                        r=r_1;
                        c=c_1;
                    else 
                        if angle_left<THRE&&angle_right>=THRE
                            back=back+1;
                            temp=zeros(s-(k+1)-1,1);
                            for q=1:s-(k+1)-1
                                right=zeros(1,3);
                                for a=1:3
                                    right(1,a)=neigh_normal((i+k-1)*n+s-k,q,a);
                                end
                                if any(right~=0)
                                    temp(q,1)=opt_func(i+k,s-k,q)+weight_1*(dist(i,i+k)+dist(i+k,i+s-1)+dist(i+s-1,i))+...
                                        weight_2*punishFun(getVecAngle(tri_normal((i-1)*n^2+(i+k-1)*n+i+s-1,:),right)*180/pi,-1,weight_2);
                                else
                                    temp(q,1)=Inf;
                                end
                            end
                            [opt_func(i,s,k),c]=min(temp);
                            
                            if opt_idx(i,k+1)==r_1
                                noback_correct=noback_correct+1;
                            else
                                if abs((opt_ref_1-opt_val_1)/opt_val_1)<=0.1
                                    noback_10=noback_10+1;
                                else
                                    if abs((opt_ref_1-opt_val_1)/opt_val_1)<=0.2
                                        noback_20=noback_20+1;
                                    else
                                        noback_more=noback_more+1;
                                    end
                                end
                            end
                            
                            if opt_idx(i+k,s-k)==c
                                back_noneed=back_noneed+1;
                            end

                            r=opt_idx(i,k+1);
                            left=zeros(1,3);
                            for a=1:3
                                left(1,a)=neigh_normal((i-1)*n+k+1,r,a);   
                            end
                            if any(left~=0)
                                opt_func(i,s,k)=opt_func(i,s,k)+opt_func(i,k+1,r)+weight_2*punishFun(getVecAngle(tri_normal((i-1)*n^2+(i+k-1)*n+i+s-1,:),left)*180/pi,-1,weight_2);
                            else
                                opt_func(i,s,k)=Inf;
                            end
                        end
                        if angle_left>=THRE&&angle_right<THRE
                            back=back+1;
                            temp=zeros((k+1)-2,1);
                            for p=1:(k+1)-2
                                left=zeros(1,3);
                                for a=1:3
                                    left(1,a)=neigh_normal((i-1)*n+k+1,p,a);   
                                end
                                if any(left~=0)
                                    temp(p,1)=opt_func(i,k+1,p)+weight_1*(dist(i,i+k)+dist(i+k,i+s-1)+dist(i+s-1,i))+...
                                        weight_2*punishFun(getVecAngle(tri_normal((i-1)*n^2+(i+k-1)*n+i+s-1,:),left)*180/pi,-1,weight_2);
                                else
                                    temp(p,1)=Inf;
                                end
                            end
                            [opt_func(i,s,k),r]=min(temp);
                            if opt_idx(i,k+1)==r
                                back_noneed=back_noneed+1;
                            end
                            if opt_idx(i+k,s-k)==c_1
                                noback_correct=noback_correct+1;
                            else
                                if abs((opt_ref_2-opt_val_2)/opt_val_2)<=0.1
                                    noback_10=noback_10+1;
                                else
                                    if abs((opt_ref_2-opt_val_2)/opt_val_2)<=0.2
                                        noback_20=noback_20+1;
                                    else
                                        noback_more=noback_more+1;
                                    end
                                end
                            end

                            c=opt_idx(i+k,s-k);
                            right=zeros(1,3);
                            for a=1:3
                                right(1,a)=neigh_normal((i+k-1)*n+s-k,c,a);   
                            end
                            if any(right~=0)
                                opt_func(i,s,k)=opt_func(i,s,k)+opt_func(i+k,s-k,c)+weight_2*punishFun(getVecAngle(tri_normal((i-1)*n^2+(i+k-1)*n+i+s-1,:),right)*180/pi,-1,weight_2);
                            else
                                opt_func(i,s,k)=Inf;
                            end
                        end
                        if angle_left<THRE&&angle_right<THRE
                            left=zeros(1,3);
                            right=zeros(1,3);
                            
                            if opt_idx(i,k+1)==r_1
                                noback_correct=noback_correct+1;
                            else
                                if abs((opt_ref_1-opt_val_1)/opt_val_1)<=0.1
                                    noback_10=noback_10+1;
                                else
                                    if abs((opt_ref_1-opt_val_1)/opt_val_1)<=0.2
                                        noback_20=noback_20+1;
                                    else
                                        noback_more=noback_more+1;
                                    end
                                end
                            end

                            if opt_idx(i+k,s-k)==c_1
                                noback_correct=noback_correct+1;
                            else
                                if abs((opt_ref_2-opt_val_2)/opt_val_2)<=0.1
                                    noback_10=noback_10+1;
                                else
                                    if abs((opt_ref_2-opt_val_2)/opt_val_2)<=0.2
                                        noback_20=noback_20+1;
                                    else
                                        noback_more=noback_more+1;
                                    end
                                end
                            end

                            p=opt_idx(i,k+1);
                            q=opt_idx(i+k,s-k);
                            for a=1:3
                                left(1,a)=neigh_normal((i-1)*n+k+1,p,a);
                                right(1,a)=neigh_normal((i+k-1)*n+s-k,q,a);
                            end
                            if any(left~=0)&&any(right~=0) 
                                opt_func(i,s,k)=opt_func(i,k+1,p)+opt_func(i+k,s-k,q)+weight_1*(dist(i,i+k)+...
                                    dist(i+k,i+s-1)+dist(i+s-1,i))+weight_2*punishFun(getVecAngle(tri_normal((i-1)*n^2+(i+k-1)*n+i+s-1,:),left)*180/pi,-1,weight_2)+...
                                    +weight_2*punishFun(getVecAngle(tri_normal((i-1)*n^2+(i+k-1)*n+i+s-1,:),right)*180/pi,-1,weight_2);
                            else
                                opt_func(i,s,k)=Inf;
                            end
                            r=p;
                            c=q;
                        end
                    end
                    face_left=face{getStoragePos(i,k+1,n),1};
                    face_right=face{getStoragePos(i+k,s-k,n),1};
                    face_temp(k,1:(k+1-2),:)=face_left(r,:,:);
                    face_temp(k,(k+1-2+1):s-2-1,:)=face_right(c,:,:);
                    face_temp(k,s-2,:)=[i,i+k,i+s-1];
                    
                    if opt_func(i,s,k)<500
                        neigh_angle((i-1)*n+s,k,1)=getVecAngle(tri_normal((i-1)*n^2+(i+k-1)*n+i+s-1,:),tri_normal((i+k-1)*n^2+(i+k+c-1)*n+i+s-1,:))*180/pi;
                        neigh_angle((i-1)*n+s,k,2)=getVecAngle(tri_normal((i-1)*n^2+(i+k-1)*n+i+s-1,:),tri_normal((i-1)*n^2+(i+r-1)*n+i+k,:))*180/pi;
                    end
                    for b=1:2
                        if neigh_angle((i-1)*n+k+1,r,b)~=0
                            count=count+1;
                            distribution(count,1)=neigh_angle((i-1)*n+k+1,r,b);
                        end
                        if neigh_angle((i+k-1)*n+s-k,c,b)~=0
                            count=count+1;
                            distribution(count,1)=neigh_angle((i+k-1)*n+s-k,c,b);
                        end
                    end
                else
                    if k==1
                        total=total+1;

                        temp=zeros(s-(k+1)-1,1);
                        right=zeros(1,3);

                        for q=1:s-(k+1)-1
                            for a=1:3
                                right(1,a)=neigh_normal((i+k-1)*n+s-k,q,a);
                            end
                            if any(right~=0)
                                temp(q,1)=opt_func(i+k,s-k,q)+...
                                    weight_2*punishFun(getVecAngle(tri_normal((i-1)*n^2+(i+k-1)*n+i+s-1,:),right)*180/pi,-1,weight_2);
                            else
                                temp(q,1)=Inf;
                                continue;
                            end
                        end
                        [opt_val,min_id]=min(temp);

                        right_opt_normal=zeros(1,3);
                        for a=1:3
                            right_opt_normal(1,a)=neigh_normal((i+k-1)*n+s-k,opt_idx(i+k,s-k),a);
                        end

                        opt_ref=opt_func(i+k,s-k,opt_idx(i+k,s-k))+weight_2*punishFun(getVecAngle(tri_normal((i-1)*n^2+(i+k-1)*n+i+s-1,:),right_opt_normal)*180/pi,-1,weight_2);

                        angle_right=getVecAngle(tri_normal((i-1)*n^2+(i+k-1)*n+i+s-1,:),right_opt_normal)*180/pi;
                        if angle_right>=THRE
                            back=back+1;
                            if min_id==opt_idx(i+k,s-k)
                                back_noneed=back_noneed+1;
                            end
                            opt_func(i,s,k)=opt_val+weight_1*(dist(i,i+k)+dist(i+k,i+s-1)+dist(i+s-1,i));
                            if any(neigh_face_normal(i,:)~=0) 
                                opt_func(i,s,k)=opt_func(i,s,k)+weight_2*punishFun(getVecAngle(tri_normal((i-1)*n^2+(i+k-1)*n+i+s-1,:),neigh_face_normal(i,:))*180/pi,is_constrained(i,1),weight_2);
                            end
                        else
                            if any(neigh_face_normal(i,:)~=0) 
                                opt_func(i,s,k)=opt_func(i,s,k)+weight_2*punishFun(getVecAngle(tri_normal((i-1)*n^2+(i+k-1)*n+i+s-1,:),neigh_face_normal(i,:))*180/pi,is_constrained(i,1),weight_2);
                            end
                            right=zeros(1,3);

                            if min_id==opt_idx(i+k,s-k)
                                noback_correct=noback_correct+1;
                            else
                                if abs((opt_ref-opt_val)/opt_val)<=0.1
                                    noback_10=noback_10+1;
                                else
                                    if abs((opt_ref-opt_val)/opt_val)<=0.2
                                        noback_20=noback_20+1;
                                    else
                                        noback_more=noback_more+1;
                                    end
                                end
                            end

                            q=opt_idx(i+k,s-k);
                            for a=1:3
                                right(1,a)=neigh_normal((i+k-1)*n+s-k,q,a);
                            end
                            if any(right~=0)
                                opt_func(i,s,k)=opt_func(i,s,k)+opt_func(i+k,s-k,q)+weight_1*(dist(i,i+k)+dist(i+k,i+s-1)+dist(i+s-1,i))+...
                                    weight_2*punishFun(getVecAngle(tri_normal((i-1)*n^2+(i+k-1)*n+i+s-1,:),right)*180/pi,-1,weight_2);
                            else
                                opt_func(i,s,k)=Inf;
                            end
                            min_id=q;
                        end
                        face_right=face{getStoragePos(i+k,s-k,n),1};
                        face_temp(k,(k+1-2+1):s-2-1,:)=face_right(min_id,:,:);
                        face_temp(k,s-2,:)=[i,i+k,i+s-1];
                        if opt_func(i,s,k)<500
                            neigh_angle((i-1)*n+s,k,1)=getVecAngle(tri_normal((i-1)*n^2+(i+k-1)*n+i+s-1,:),tri_normal((i+k-1)*n^2+(i+k+min_id-1)*n+i+s-1,:))*180/pi;
                        end
                        for b=1:2
                            if neigh_angle((i+k-1)*n+s-k,min_id,b)~=0
                                count=count+1;
                                distribution(count,1)=neigh_angle((i+k-1)*n+s-k,min_id,b);
                            end
                        end
                    else
                        total=total+1;
                        temp=zeros((k+1)-2,1);
                        left=zeros(1,3);

                        for p=1:(k+1)-2
                            for a=1:3
                                left(1,a)=neigh_normal((i-1)*n+k+1,p,a);   
                            end
                            if any(left~=0)
                                temp(p,1)=opt_func(i,k+1,p)+...
                                    weight_2*punishFun(getVecAngle(tri_normal((i-1)*n^2+(i+k-1)*n+i+s-1,:),left)*180/pi,-1,weight_2);
                            else
                                temp(p,1)=Inf;
                                continue;
                            end
                        end
                        [opt_val,min_id]=min(temp);

                        left_opt_normal=zeros(1,3);
                        for a=1:3
                            left_opt_normal(1,a)=neigh_normal((i-1)*n+k+1,opt_idx(i,k+1),a);
                        end

                        opt_ref=opt_func(i,k+1,opt_idx(i,k+1))+weight_2*punishFun(getVecAngle(tri_normal((i-1)*n^2+(i+k-1)*n+i+s-1,:),left_opt_normal)*180/pi,-1,weight_2);

                        angle_left=getVecAngle(tri_normal((i-1)*n^2+(i+k-1)*n+i+s-1,:),left_opt_normal)*180/pi;
                        if angle_left>=THRE
                            back=back+1;
                            if min_id==opt_idx(i,k+1)
                                back_noneed=back_noneed+1;
                            end
                            opt_func(i,s,k)=opt_val+weight_1*(dist(i,i+k)+dist(i+k,i+s-1)+dist(i+s-1,i));
                            if any(neigh_face_normal(i+s-2,:)~=0) 
                                opt_func(i,s,k)=opt_func(i,s,k)+weight_2*punishFun(getVecAngle(tri_normal((i-1)*n^2+(i+k-1)*n+i+s-1,:),neigh_face_normal(i+s-2,:))*180/pi,is_constrained(i+s-2,1),weight_2);
                            end
                        else
                            if any(neigh_face_normal(i+s-2,:)~=0) 
                                opt_func(i,s,k)=opt_func(i,s,k)+weight_2*punishFun(getVecAngle(tri_normal((i-1)*n^2+(i+k-1)*n+i+s-1,:),neigh_face_normal(i+s-2,:))*180/pi,is_constrained(i+s-2,1),weight_2);
                            end
                            left=zeros(1,3);
                            
                            if min_id==opt_idx(i,k+1)
                                noback_correct=noback_correct+1;
                            else
                                if abs((opt_ref-opt_val)/opt_val)<=0.1
                                    noback_10=noback_10+1;
                                else
                                    if abs((opt_ref-opt_val)/opt_val)<=0.2
                                        noback_20=noback_20+1;
                                    else
                                        noback_more=noback_more+1;
                                    end
                                end
                            end

                            p=opt_idx(i,k+1);
                            for a=1:3
                                left(1,a)=neigh_normal((i-1)*n+k+1,p,a);   
                            end
                            if any(left~=0)
                                opt_func(i,s,k)=opt_func(i,s,k)+opt_func(i,k+1,p)+weight_1*(dist(i,i+k)+dist(i+k,i+s-1)+dist(i+s-1,i))+...
                                    weight_2*punishFun(getVecAngle(tri_normal((i-1)*n^2+(i+k-1)*n+i+s-1,:),left)*180/pi,-1,weight_2);
                            else
                                opt_func(i,s,k)=Inf;
                            end
                            min_id=p;
                        end
                        face_left=face{getStoragePos(i,k+1,n),1};
                        face_temp(k,1:(k+1-2),:)=face_left(min_id,:,:);
                        face_temp(k,s-2,:)=[i,i+k,i+s-1];
                        if opt_func(i,s,k)<500
                            neigh_angle((i-1)*n+s,k,1)=getVecAngle(tri_normal((i-1)*n^2+(i+k-1)*n+i+s-1,:),tri_normal((i-1)*n^2+(i+min_id-1)*n+i+k,:))*180/pi;
                        end
                        for b=1:2
                            if neigh_angle((i-1)*n+k+1,min_id,b)~=0
                                count=count+1;
                                distribution(count,1)=neigh_angle((i-1)*n+k+1,min_id,b);
                            end
                        end
                    end
                end
                face{getStoragePos(i,s,n),1}=face_temp;
                neigh_normal((i-1)*n+s,k,:)=tri_normal((i-1)*n^2+(i+k-1)*n+i+s-1,:);
                if s==n
                    if any(neigh_face_normal(n,:)~=0)
                        opt_func(i,s,k)=opt_func(i,s,k)+weight_2*punishFun(getVecAngle(tri_normal((i-1)*n^2+(i+k-1)*n+i+s-1,:),neigh_face_normal(n,:))*180/pi,is_constrained(n,1),weight_2);
                    end
                end
            end
            [~,opt_idx(i,s)]=min(opt_func(i,s,1:s-2));
        end
    end

    opt_val=opt_func(1,n,opt_idx(1,n));
    opt_face=zeros(n-2,3);
    opt_face_cell=face{getStoragePos(1,n,n),1};
    for i=1:n-2
        opt_face(i,:)=opt_face_cell(opt_idx(1,n),i,:);
    end
end

function [threshold]=getThreshold(distribution,p)
    distribution=sort(distribution,'descend');
    data_num=ceil(p*size(distribution,1));
    threshold=distribution(data_num,1);
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

function [pos]=getStoragePos(i,s,n)
    pos=(s-3)*((n-2)+(n-(s-2)))/2+i;
end
