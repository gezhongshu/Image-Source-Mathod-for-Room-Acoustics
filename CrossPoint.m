function [point] = CrossPoint(imagepoint,rcv,plane,g)
% This function computes the crossing point between a line and a plane
k = imagepoint - rcv;
b = rcv;
%track = [k(2)+k(3),-k(1),-k(1),-(k(2)*b(1)-k(1)*b(2)+k(3)*b(1)-k(1)*b(3))]
t = -(dot(plane(g,1:3),b(1:3))+plane(g,4))/dot(plane(g,1:3),k(1:3));
point = [k(1)*t+b(1),k(2)*t+b(2),k(3)*t+b(3)];