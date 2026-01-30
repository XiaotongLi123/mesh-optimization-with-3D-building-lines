function []=writeLineOBJ(filename,line_point,line_sign)
    % Write the selected line elements (line_sign == 1) to file
    line_num=sum(line_sign);
    file=fopen(filename,'w');
    for i=1:size(line_point,1)/2
        if line_sign(i,1)==1
            fprintf(file,'v %f %f %f\n',line_point(2*i-1,1),line_point(2*i-1,2),line_point(2*i-1,3));
            fprintf(file,'v %f %f %f\n',line_point(2*i,1),line_point(2*i,2),line_point(2*i,3));
        end
    end
    for i=1:line_num
        fprintf(file,'l %d %d\n',2*i-1,2*i);
    end
    fclose(file);
end