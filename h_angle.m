function [coef] = h_angle(coef,plane,vec1,vec2,vec3,n)

v1 = dot(plane(n,1:3),vec1/norm(vec1))*vec1;
v2 = dot(plane(n,1:3),vec2/norm(vec2))*vec2;
v = v1 + v2;

if v(1) == 0 && v(2) == 0
    coef(n,2) = 0;
else
    temp = dot(plane(n,1:2),[1,0])/norm(plane(n,1:2));
    temp1 = cross(v,vec1);
    if temp1(1)<0 || temp1(2)<0 || temp1(3)<0
        coef(n,2) = acos(temp);
    else
        coef(n,2) = -acos(temp);
    end
    
end

% if plane(n,1:2) == 0
%     coef(n,2) = 0;
% else
%     temp = dot(plane(n,1:2),[1,0])/norm(plane(n,1:2));
%     if plane(n,2) >= 0
%         coef(n,2) = acos(temp);
%     else
%         coef(n,2) = -acos(temp);
%     end
% end