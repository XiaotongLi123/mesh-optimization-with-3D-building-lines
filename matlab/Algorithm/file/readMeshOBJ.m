function [point,face]=readMeshOBJ(filename)
    point_ori=zeros(9999999,3);
    point_num=0;
    face_ori=zeros(9999999,3);
    face_num=0;
    file=fopen(filename,'r');
    while ~feof(file)
        str=fgetl(file);

        s=regexp(str,' ','split');
        sign=char(s{1});

        if sign=='#'
            continue;
        end

        if sign=='v'
            point_num=point_num+1;
            x=str2num(char(s{2}));
            y=str2num(char(s{3}));
            z=str2num(char(s{4}));
            point_ori(point_num,:)=[x,y,z];
        end

        if sign=='f'
            face_type=size(s,2)-1;
            if face_type<3
                continue;
            else 
                if face_type==3
                    face_num=face_num+1;
                    v1=str2num(char(s{2}));
                    v2=str2num(char(s{3}));
                    v3=str2num(char(s{4}));
                    face_ori(face_num,:)=[v1,v2,v3];
                else
                    for i=2:face_type-1
                        face_num=face_num+1;
                        v1=str2num(char(s{2}));
                        v2=str2num(char(s{i}));
                        v3=str2num(char(s{i+1}));
                        face_ori(face_num,:)=[v1,v2,v3];
                    end
                end
            end
            
        end
    end
    fclose(file);

    point=point_ori(1:point_num,:);
    face=face_ori(1:face_num,:);
end


% function [point,face]=readMeshOBJ(filename)
%     point_ori=zeros(999999,3);
%     point_num=0;
%     face_ori=zeros(999999,3);
%     face_num=0;
%     file=fopen(filename,'r');
%     data=textscan(file,'%s%f%f%f');
%     fclose(file);
% 
%     sign=data{1,1};
%     col_2=data{1,2};
%     col_3=data{1,3};
%     col_4=data{1,4};
%     for i=1:size(sign,1)
%         if sign{i,1}(1,1)=='v'
%             point_num=point_num+1;
%             point_ori(point_num,:)=[col_2(i,1),col_3(i,1),col_4(i,1)];
%         end
%         if sign{i,1}(1,1)=='f'
%             face_num=face_num+1;
%             face_ori(face_num,:)=[col_2(i,1),col_3(i,1),col_4(i,1)];
%         end
%     end
%     point=point_ori(1:point_num,:);
%     face=face_ori(1:face_num,:);
% end