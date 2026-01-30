function [node]=getNode(line,interval)
    constrained_line_length=point2pointDist(line(1,:),line(2,:));
    segment_num=round(constrained_line_length/interval);

    if segment_num==0
        segment_num=1;
    end

    interval=point2pointDist(line(1,:),line(2,:))/segment_num;
    node_num=segment_num+1;

    [a,b,c]=LineEquation(line(1,:),line(2,:));

    node=zeros(node_num,3);
    for i=1:node_num-1
        t=(i-1)*interval/sqrt(a^2+b^2+c^2);
        node(i,1)=line(1,1)+a*t;
        node(i,2)=line(1,2)+b*t;
        node(i,3)=line(1,3)+c*t;
    end
    node(end,:)=line(2,:);
end