function [angle]=getVecAngle(vec_1,vec_2)
    % MATLAB may return complex values when vector angles are too small or too large
    % The computed angles are given in radians
    angle=real(acos(dot(vec_1,vec_2)/(norm(vec_1)*norm(vec_2))));
end