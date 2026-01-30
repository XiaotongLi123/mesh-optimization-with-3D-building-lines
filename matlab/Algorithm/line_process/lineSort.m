function [line_sort]=lineSort(line)
    line_num=size(line,1)/2;
    length=zeros(line_num,1);
    for i=1:line_num
        length(i,1)=point2pointDist(line(2*i-1,:),line(2*i,:));
    end
    [~,sort_ind]=sort(length,'descend');
    line_sort=zeros(line_num*2,3);

    for i=1:line_num
        line_sort(2*i-1:2*i,:)=line(2*sort_ind(i,1)-1:2*sort_ind(i,1),:);
    end
end