function [coef] = h_angle (coef,plane,n)
if plane(n,1:2) == 0
    coef(n,2) = 0;
else
    temp = dot(plane(n,1:2),[1,0])/norm(plane(n,1:2));
    if plane(n,2) >= 0
        coef(n,2) = acos(temp);
    else
        coef(n,2) = -acos(temp);
    end
end