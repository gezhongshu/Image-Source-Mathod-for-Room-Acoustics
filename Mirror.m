function [image]  = Mirror(wall,vertex,src,n,plane)
v = src - vertex(wall(n,2),:);
v1 = vertex(wall(n,1),:) - vertex(wall(n,2),:);
v2 = vertex(wall(n,3),:) - vertex(wall(n,2),:);
v1 = v1/norm(v1);
v2 = v2/norm(v2);
%-----------------------------------------------------------------------%
% if two vectors are not orthogornal, create a vector in the plane that can
% form orthogornal to the other vertor This part should be rechecked
% carefully!!!
if dot(v1,v2)~=0
    %     error('Unorthogornal');
    count1 = plane(n,1:3);
    count1 = count1(count1~=0);
    count2 = v1(v1~=0);
    if size(count1,2) == 1 && size(count2,2) == 1
        pos1 = find(plane(n,1:3)~=0);
        pos2 = find(v1~=0);
        v2 = ones(1,3);
        v2(pos1) = 0;
        v2(pos2) = 0;
    else
        temp = [plane(n,2:3);v1(2:3)]\[-plane(n,1);-v1(1)];
        temp = temp';
        v2 = [1 temp];
        v2 = v2/norm(v2);
    end
end
%-----------------------------------------------------------------------%
temp1 = dot(v,v1)*v1;
temp2 = dot(v,v2)*v2;
image=(-(v-(temp1+temp2))+temp1+temp2)+vertex(wall(n,2),:);


% 
% v3 = temp1+temp2;
% v4 = v - v3;
% dot(v3,v4)