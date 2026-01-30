function []=writePointOBJ(filename,point)
    file=fopen(filename,'w');
    for i=1:size(point,1)
        fprintf(file,'v %f %f %f\n',point(i,1),point(i,2),point(i,3));
    end
    fclose(file);
end