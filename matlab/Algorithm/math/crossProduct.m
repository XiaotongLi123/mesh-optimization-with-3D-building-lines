function [result_vector]=crossProduct(vector)
    result_vector=zeros(1,3);
    result_vector(1,1)=vector(1,2)*vector(2,3)-vector(2,2)*vector(1,3);
    result_vector(1,2)=vector(1,3)*vector(2,1)-vector(2,3)*vector(1,1);
    result_vector(1,3)=vector(1,1)*vector(2,2)-vector(2,1)*vector(1,2);
end