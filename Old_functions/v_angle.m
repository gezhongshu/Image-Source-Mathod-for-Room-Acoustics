function [coef] = v_angle(coef,plane,vec1,vec2,vec3,n)

temp = dot(plane(n,1:3),vec3)/(norm(plane(n,1:3))*norm(vec3));
v1 = dot(plane(n,1:3),vec1/norm(vec1))*vec1;
v2 = dot(plane(n,1:3),vec3/norm(vec3))*vec3;
v = v1 + v2;

temp1 = dot(v,vec1);
temp2 = cross(v,vec1);

if temp2(1)==0 && temp2(2)==0 && temp2(3)==0
    coef(n,1) = 0;
else
    if temp1<0
        if temp2(1)>0 || temp2(2)>0 || temp2(3)>0
            coef(n,1) = -pi-(pi/2 - acos(temp));
        else
            coef(n,1) = pi-(pi/2 - acos(temp));
        end
    else
        coef(n,1) = pi/2 - acos(temp);
    end
end

% if plane(n,1)<0 && plane(n,3)>=0
%     coef(n,1) = pi-(pi/2 - acos(temp));
% else if plane(n,1)<0 && plane(n,3)<0
%         coef(n,1) = -pi-(pi/2 - acos(temp));
%     end
% end
% if plane(n,1)>=0
%     coef(n,1) = pi/2 - acos(temp);
% end