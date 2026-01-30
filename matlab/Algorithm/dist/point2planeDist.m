function [dist]=point2planeDist(point,plane)
    %plane为1-by-4矩阵，记录平面方程的参数
    dist=abs(plane(1,1)*point(1,1)+plane(1,2)*point(1,2)+plane(1,3)*point(1,3)+plane(1,4))/sqrt(plane(1,1)*plane(1,1)+plane(1,2)*plane(1,2)+plane(1,3)*plane(1,3));
end