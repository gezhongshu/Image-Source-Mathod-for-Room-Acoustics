function [image]  = Mirror(wall,vertex,src,n,plane)
v = src - vertex(wall(n,2),:);
v1 = vertex(wall(n,1),:) - vertex(wall(n,2),:);
v2 = vertex(wall(n,3),:) - vertex(wall(n,2),:);
v1 = v1/norm(v1);
v2 = v2/norm(v2);
%-----------------------------------------------------------------------%
% if two vectors are not orthogornal, create a vector in the plane that can
% form orthogornal to the other vertor
if dot(v1,v2)~=0
v2 = cross(plane(n,1:3),v1);
v2 = v2/norm(v2);
end
%-----------------------------------------------------------------------%
temp1 = dot(v,v1)*v1;
temp2 = dot(v,v2)*v2;
image=(-(v-(temp1+temp2))+temp1+temp2)+vertex(wall(n,2),:);


% 
% v3 = temp1+temp2;
% v4 = v - v3;
% dot(v3,v4)
