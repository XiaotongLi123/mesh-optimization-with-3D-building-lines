function [a, b, c] = LineEquation(V1, V2)
    % V1 and V2 are 1-by-3 vectors storing two points on a line
    % Compute the parametric equation of the line passing through V1 and V2
    
    x1 = V1(1,1);
    x2 = V2(1,1);
    y1 = V1(1,2);
    y2 = V2(1,2);
    z1 = V1(1,3);
    z2 = V2(1,3);

    % Compute line direction vector (a, b, c) for the parametric equation:
    % x = x1 + a * t
    % y = y1 + b * t
    % z = z1 + c * t
    if x1 == x2 && y1 == y2 && z1 == z2
        a = 0;
        b = 0;
        c = 0; % Two points are identical, direction vector is zero
    else
        norm_val = sqrt((x2 - x1)^2 + (y2 - y1)^2 + (z2 - z1)^2);
        a = (x2 - x1) / norm_val;
        b = (y2 - y1) / norm_val;
        c = (z2 - z1) / norm_val;
    end
end