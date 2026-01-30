function [is_exist_point,exist_point_id]=isExistPoint(point_set,this_point)
    pts_num=size(point_set,1);
    is_exist_point=0;
    exist_point_id=zeros(pts_num,1);
    count=0;
    for i=1:pts_num
        if point2pointDist(point_set(i,:),this_point)<0.001
            is_exist_point=1;
            count=count+1;
            exist_point_id(count,1)=i;
        end
    end
    exist_point_id=exist_point_id(1:count,1);
end