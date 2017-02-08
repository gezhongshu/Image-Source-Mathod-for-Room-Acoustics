function [flag] = BoundaryJudge(wall,vertex,point,n)
% test if the collision point is in the polygon of walls. if flag = 1, the
% point is in the polygon. If flag = 0, the point is not in the polygon.

comp = [];
flag = 1;
th = zeros(1,size(wall,2));
endpoint = 0;
for m = 1:1:size(wall,2)
    % ensure the vertices number is not zero
    % test if it is the last point to do the last product.
    if (m == size(wall,2) && wall(n,m)~=0)||(wall(n,m) ~= 0 && wall(n,m+1)==0)
        v1 = vertex(wall(n,m),:)-point(1:3);
        v2 = vertex(wall(n,1),:)-point(1:3);
        temp = cross(v1,v2);
        theta = acos(dot(v1,v2)/(norm(v1)*norm(v2)));
        endpoint = 1;
    else if wall(n,m) ~= 0 && wall(n,m+1) ~= 0;
            v1 = vertex(wall(n,m),:)-point(1:3);
            v2 = vertex(wall(n,m+1),:)-point(1:3);
            temp = cross(v1,v2);
            theta = acos(dot(v1,v2)/(norm(v1)*norm(v2)));
        end
    end
    % check if comparing vector is empty, if it is not empty, do the
    % judgement
    if isempty(comp)
        comp = temp;
    end
    % decide if two vectors are in the same direction
    if dot(temp,comp)<0
        th(m) = -theta;
    else
        th(m) = theta;
    end
    if endpoint == 1;
        break
    end
end
if abs(sum(th)) <= 1.0000e-04
    flag = 0;
end