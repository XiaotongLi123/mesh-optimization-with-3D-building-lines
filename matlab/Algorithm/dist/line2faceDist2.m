function [l2f_dist]=line2faceDist2(tri,line)
    % Perpendicular distance from a line to a triangular face
    plane=getPlanePara(tri);

    p2f_dist_1=point2planeDist(line(1,:),plane);
    p2f_dist_2=point2planeDist(line(2,:),plane);

    l2f_dist=max(p2f_dist_1,p2f_dist_2);
end

function [para]=getPlanePara(tri)
    V1=tri(1,:);
    V2=tri(2,:);
    V3=tri(3,:);
            
    % Plane equation
    %|x-x1    y-y1    z-z1 |
    %|x2-x1   y2-y1   z2-z1| = 0
    %|x3-x1   y3-y1   z3-z1|
    
    %s11x+s12y+s13z+(-s11x1-s12y1-s13z1)=0
    a=det([V2(2)-V1(2),V2(3)-V1(3);V3(2)-V1(2),V3(3)-V1(3)]);
    b=-det([V2(1)-V1(1),V2(3)-V1(3);V3(1)-V1(1),V3(3)-V1(3)]);
    c=det([V2(1)-V1(1),V2(2)-V1(2);V3(1)-V1(1),V3(2)-V1(2)]);
    d=-a*V1(1)-b*V1(2)-c*V1(3);

    para=[a,b,c,d];
end