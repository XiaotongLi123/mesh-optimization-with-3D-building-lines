function [point,line]=readLineOBJ(filename)
    point_ori=zeros(999999,3);
    point_num=0;
    line_ori=zeros(999999,2);
    line_num=0;
    file=fopen(filename,'r');
    data=textscan(file,'%s%f%f%f');
    fclose(file);

    sign=data{1,1};
    col_2=data{1,2};
    col_3=data{1,3};
    col_4=data{1,4};
    for i=1:size(sign,1)
        if sign{i,1}(1,1)=='v'
            point_num=point_num+1;
            point_ori(point_num,:)=[col_2(i,1),col_3(i,1),col_4(i,1)];
        end
        if sign{i,1}(1,1)=='l'
            line_num=line_num+1;
            line_ori(line_num,:)=[col_2(i,1),col_3(i,1)];
        end
    end
    point=point_ori(1:point_num,:);
    line=line_ori(1:line_num,:);
end