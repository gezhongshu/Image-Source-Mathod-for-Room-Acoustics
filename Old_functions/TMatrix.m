function [m] = TMatrix(plane,wnum,vec1,vec2,vec3,Opoint,n)
% coef = zeros(wnum,3); % coefficients container
% matrix = zeros(4,4,wnum);
% for n = 1:1:wnum
%     % (1)Calculate the angle of vector and the xy plane
%     coef = v_angle(coef,plane,vec,n);
%     % cos(coef(5,1)) % 6.1232e-17
%     % (2)Calculate the angle between the plane vector projection on the xy
%     % plane and the x axis
%     coef = h_angle (coef,plane,n);
%     % (3)Calculate the distance from the original point to the plane
%     coef(n,3) = abs(dot(plane(n,1:3),Opoint)+plane(n,4))/norm(plane(n,1:3)-Opoint);
%     % Compute the transforming matrix
%     s2 = sin(coef(n,2));
%     c2 = cos(coef(n,2));
%     s1 = sin(coef(n,1));
%     c1 = cos(coef(n,1));
%     s21 = sin(2*coef(n,1));
%     c21 = cos(2*coef(n,1));
%     m = [s2^2-c2^2*c21,     -c2*s2*(c21+1),     -c2*s21,        0; ...
%         -c2*s2*(c21+1),     -s2^2*c21+c2^2,     -s2*s21,        0; ...
%         -c2*s21,            -s2*s21,            c21,            0; ...
%         2*coef(n,3)*c1*c2,  2*coef(n,3)*c1*s2,  2*coef(n,3)*s1  1];
%     matrix(:,:,n) = m;
% end

coef = zeros(wnum,3); % coefficients container
% (1)Calculate the angle of vector and the xy plane
coef = v_angle(coef,plane,vec1,vec2,vec3,n);
% cos(coef(5,1)) % 6.1232e-17
% (2)Calculate the angle between the plane vector projection on the xy
% plane and the x axis
coef = h_angle(coef,plane,vec1,vec2,vec3,n);
% (3)Calculate the distance from the original point to the plane
coef(n,3) = abs(dot(plane(n,1:3),Opoint)+plane(n,4))/norm(plane(n,1:3));
% Compute the transforming matrix
s2 = sin(coef(n,2));
c2 = cos(coef(n,2));
s1 = sin(coef(n,1));
c1 = cos(coef(n,1));
s21 = sin(2*coef(n,1));
c21 = cos(2*coef(n,1));
m = [s2^2-c2^2*c21,     -c2*s2*(c21+1),     -c2*s21,        0; ...
    -c2*s2*(c21+1),     -s2^2*c21+c2^2,     -s2*s21,        0; ...
    -c2*s21,            -s2*s21,            c21,            0; ...
    2*coef(n,3)*c1*c2,  2*coef(n,3)*c1*s2,  2*coef(n,3)*s1  1];