function [dist]=point2pointDist(point_1,point_2)
    dist=sqrt((point_2(1,1)-point_1(1,1))^2+(point_2(1,2)-point_1(1,2))^2+(point_2(1,3)-point_1(1,3))^2);
end