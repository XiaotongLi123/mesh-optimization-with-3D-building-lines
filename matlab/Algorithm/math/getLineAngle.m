function [angle]=getLineAngle(vec_1,vec_2)
    % MATLAB may return complex values when vector angles are too small or too large
    % Computed angles are expressed in radians
    angle=real(acos(dot(vec_1,vec_2)/(norm(vec_1)*norm(vec_2))));
    if angle>pi/2
        angle=pi-angle;
    end
end