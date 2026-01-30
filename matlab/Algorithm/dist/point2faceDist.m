function [p2f_dist,nearest_point]=point2faceDist(point,triangle)
    % Distance from a point to a triangular face
    V1=triangle(1,:);
    V2=triangle(2,:);
    V3=triangle(3,:);
            
    % Plane equation
    %|x-x1    y-y1    z-z1 |
    %|x2-x1   y2-y1   z2-z1| = 0
    %|x3-x1   y3-y1   z3-z1|
    
    %s11x+s12y+s13z+(-s11x1-s12y1-s13z1)=0
    s11=det([V2(2)-V1(2),V2(3)-V1(3);V3(2)-V1(2),V3(3)-V1(3)]);
    s12=-det([V2(1)-V1(1),V2(3)-V1(3);V3(1)-V1(1),V3(3)-V1(3)]);
    s13=det([V2(1)-V1(1),V2(2)-V1(2);V3(1)-V1(1),V3(2)-V1(2)]);
    a_t=s11;
    b_t=s12;
    c_t=s13;

    % Equation of the perpendicular line
    %x=point(1,1)+a_t*t
    %y=point(1,2)+b_t*t
    %z=point(1,3)+c_t*t

    t0=-((point(1,1)-V1(1))*s11+(point(1,2)-V1(2))*s12+(point(1,3)-V1(3))*s13)/(a_t*s11+b_t*s12+c_t*s13);
    x0=point(1,1)+a_t*t0;
    y0=point(1,2)+b_t*t0;
    z0=point(1,3)+c_t*t0;

    if isInTri(triangle,[x0,y0,z0]) 
        p2f_dist=abs(s11*point(1,1)+s12*point(1,2)+s13*point(1,3)-s11*V1(1)-s12*V1(2)...
            -s13*V1(3))/sqrt(s11^2+s12^2+s13^2);
        nearest_point=[x0,y0,z0];
    else
        nearest_point_candidate=zeros(3,3);
        [dist_1,nearest_point_candidate(1,:)]=point2segmentDist(point,[V2;V3]);
        [dist_2,nearest_point_candidate(2,:)]=point2segmentDist(point,[V3;V1]);
        [dist_3,nearest_point_candidate(3,:)]=point2segmentDist(point,[V1;V2]);
        [p2f_dist,nearest_id]=min([dist_1;dist_2;dist_3]);
        nearest_point=nearest_point_candidate(nearest_id,:);
    end
end