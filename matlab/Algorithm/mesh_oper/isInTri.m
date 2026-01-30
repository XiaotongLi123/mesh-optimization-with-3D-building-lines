function [in]=isInTri(triangle,point)
    %判断三维空间点是否在三角平面上
    V1=triangle(1,:);
    V2=triangle(2,:);
    V3=triangle(3,:);
    x0=point(1,1);
    y0=point(1,2);
    z0=point(1,3);

    %计算三角形三条直线的方程
    [a1,b1,c1]=LineEquation(V2,V3);
    [a2,b2,c2]=LineEquation(V1,V3);
    [a3,b3,c3]=LineEquation(V1,V2);
    %计算点到直线距离
    D=sqrt((b1*(V1(3)-V2(3))-c1*(V1(2)-V2(2)))^2+(a1*(V1(3)-V2(3))-c1*(V1(1)-V2(1)))^2+...
        (a1*(V1(2)-V2(2))-b1*(V1(1)-V2(1)))^2)/sqrt(a1^2+b1^2+c1^2);
    D1=sqrt((b1*(z0-V2(3))-c1*(y0-V2(2)))^2+(a1*(z0-V2(3))-c1*(x0-V2(1)))^2+...
        (a1*(y0-V2(2))-b1*(x0-V2(1)))^2)/sqrt(a1^2+b1^2+c1^2);
    D2=sqrt((b2*(z0-V3(3))-c2*(y0-V3(2)))^2+(a2*(z0-V3(3))-c2*(x0-V3(1)))^2+...
        (a2*(y0-V3(2))-b2*(x0-V3(1)))^2)/sqrt(a2^2+b2^2+c2^2);
    D3=sqrt((b3*(z0-V1(3))-c3*(y0-V1(2)))^2+(a3*(z0-V1(3))-c3*(x0-V1(1)))^2+...
        (a3*(y0-V1(2))-b3*(x0-V1(1)))^2)/sqrt(a3^2+b3^2+c3^2);

    %计算三条线段的长度
    L1=sqrt((V2(1)-V3(1))^2+(V2(2)-V3(2))^2+(V2(3)-V3(3))^2);
    L2=sqrt((V1(1)-V3(1))^2+(V1(2)-V3(2))^2+(V1(3)-V3(3))^2);
    L3=sqrt((V2(1)-V1(1))^2+(V2(2)-V1(2))^2+(V2(3)-V1(3))^2);
    Area=D*L1/2;
    Area_1=D1*L1/2;
    Area_2=D2*L2/2;
    Area_3=D3*L3/2;
    error=Area-(Area_1+Area_2+Area_3);
    
    %考虑舍入误差，error并不严格等于0
    if abs(error)<Area*0.0000001
        in=1;
    else
        in=0;
    end
end