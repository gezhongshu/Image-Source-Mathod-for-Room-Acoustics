function [coef] = v_angle (coef,plane,vec,n)
temp = dot(plane(n,1:3),vec)/norm(plane(n,1:3));
if plane(n,1)<0 && plane(n,3)>=0
    coef(n,1) = pi-(pi/2 - acos(temp));
else if plane(n,1)<0 && plane(n,3)<0
        coef(n,1) = -pi-(pi/2 - acos(temp));
    end
end
if plane(n,1)>=0
    coef(n,1) = pi/2 - acos(temp);
end