function [plane] = TPlane(wall,wnum,vertex,src)
% transform data form of wall plus vertex to corresponding wall functions
% in 3D domain. 
% wall:     arrays including vertexes each wall includes
% vertex:   arrays including each vertex's value
% src:      source point that is going to be mirrored
plane = zeros(wnum,4);
for n = 1:1:wnum
    p1 = vertex(wall(n,2),:) - vertex(wall(n,1),:);
    p2 = vertex(wall(n,3),:) - vertex(wall(n,1),:);
    A = p1(2)*p2(3)-p2(2)*p1(3);
    B = p1(3)*p2(1)-p2(3)*p1(1);
    C = p1(1)*p2(2)-p2(1)*p1(2);
    D = -(A*vertex(wall(n,1),1)+B*vertex(wall(n,1),2)+C*vertex(wall(n,1),3));
    plane(n,:) = [A B C D];
    % here we have to test the directions of the plane vectors and set them
    % go backwards the source point. This is a must for later angle
    % calculation in the axis transformation part
    testvect = vertex(wall(n,1),:)-src(1:3);
    temp = dot(testvect,plane(n,1:3));
    if temp < 0
        plane(n,1:3) = -plane(n,1:3);
    end
end