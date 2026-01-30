function [area]=getTriArea(triangle)
    a=point2pointDist(triangle(2,:),triangle(3,:));
    b=point2pointDist(triangle(1,:),triangle(3,:));
    c=point2pointDist(triangle(1,:),triangle(2,:));
    p=(a+b+c)/2;
    area=real(sqrt(p*(p-a)*(p-b)*(p-c)));
end