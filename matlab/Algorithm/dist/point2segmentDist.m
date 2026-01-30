function [dist,nearest_point]=point2segmentDist(point,line)
    [a,b,c]=LineEquation(line(1,:),line(2,:));
    d=-(a*point(1,1)+b*point(1,2)+c*point(1,3));
    t=-(a*line(1,1)+b*line(1,2)+c*line(1,3)+d)/(a^2+b^2+c^2);
    x=line(1,1)+a*t;
    y=line(1,2)+b*t;
    z=line(1,3)+c*t;
    vector_1=point-line(1,:);
    vector_2=[x,y,z]-point;
    vector_3=vector_1+vector_2;
    vector_4=line(2,:)-line(1,:);
    vector_5=line(2,:)-point;
    for i=1:3
        if vector_4(1,i)~=0
            m=i;
            break;
        end
    end
    lambda=vector_3(1,m)/vector_4(1,m); % Avoid NaN values due to division by zero when a vector component is zero
    if lambda>=0&&lambda<=1
        dist=norm(vector_2);
        nearest_point=[x,y,z];
    end
    if lambda<0
        dist=norm(vector_1);
        nearest_point=line(1,:);
    end
    if lambda>1
        dist=norm(vector_5);
        nearest_point=line(2,:);
    end
end