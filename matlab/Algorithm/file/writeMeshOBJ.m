function []=writeMeshOBJ(filename,point,face)
    point_num=size(point,1);
    face_num=size(face,1);
    file=fopen(filename,'w');
    for i=1:point_num
        fprintf(file,'v %f %f %f\n',point(i,1),point(i,2),point(i,3));
    end
    for i=1:face_num
        fprintf(file,'f %d %d %d\n',face(i,1),face(i,2),face(i,3));
    end
    fclose(file);
end