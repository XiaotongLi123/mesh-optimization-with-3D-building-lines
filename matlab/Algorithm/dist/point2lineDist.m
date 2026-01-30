function [dist]=point2lineDist(point,line)
    [a,b,c]=LineEquation(line(1,:),line(2,:));
    x0=point(1,1); y0=point(1,2); z0=point(1,3);
    dist=sqrt((b*(z0-line(1,3))-c*(y0-line(1,2)))^2+(a*(z0-line(1,3))-c*(x0-line(1,1)))^2+...
        (a*(y0-line(1,2))-b*(x0-line(1,1)))^2)/sqrt(a^2+b^2+c^2);
end
