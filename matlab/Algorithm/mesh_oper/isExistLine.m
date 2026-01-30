function [is_exist_line,exist_line_id]=isExistLine(line_set,line)
    % Check whether a given edge already exists in the line list
    line_num=size(line_set,1);
    is_exist_line=0;
    exist_line_id=zeros(line_num,1);
    count=0;

    for i=1:line_num
        if any(ismember(line_set(i,:),line(1,1)))&&any(ismember(line_set(i,:),line(1,2)))
            is_exist_line=1;
            count=count+1;
            exist_line_id(count,1)=i;
        end
    end
    exist_line_id=exist_line_id(1:count,1);
end