%% Initialization
clear;
Initial;

%% Generate the transforming matrix for each plane
coef = zeros(wnum,3); % coefficients container
vec = [0 0 1]; % vector referrence plane xy
matrix = zeros(4,4,wnum);
for n = 1:1:wnum
    % (1)Calculate the angle of vector and the xy plane
    coef = v_angle(coef,plane,vec,n);
    % cos(coef(5,1)) % 6.1232e-17
    % (2)Calculate the angle between the plane vector projection on the xy
    % plane and the x axis
    coef = h_angle (coef,plane,n);
    % (3)Calculate the distance from the original point to the plane
    coef(n,3) = abs(plane(n,4))/norm(plane(n,1:3));
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
    matrix(:,:,n) = m;
end

%% Mirroring processing
% Check the vailidity of the path
% Wall absorbing coefficients
% Air absorbing coefficients
N = 2;
image = cell(1,N);

for n = 1:1:N
    if n == 1
        for g = 1:1:wnum
            image{1,n}{1,g} = src * matrix(:,:,g);
            % Testing of the validity of the point begins
            for l = 1:1:wnum
                % calculate the cross point between the trajectory and the
                % relfecting walls
                point = CrossPoint(image{1,n}{1,g},rcv,plane,l);
                if l ~= g 
                    % test if there exist obstructions in the path
                    % (1)Test if the point is in the path of the two points
                    flag = InPath(point,image{1,n}{1,g},rcv);
                    if flag == 0
                        continue
                    else
                        % (2)Test if the collision point is in the polygon
                        % of walls. If the point IS in the polygon, there
                        % exist obstructions in the path and the point
                        % should be removed
                        flag = BoundaryJudge(wall,vertex,point,l);
                        if flag == 1
                            image{1,n}{1,g} = rcv;
                        end
                    end
                end
                if l == g
                    % (1)Test if the point is in the path of the two points
                    flag = InPath(point,image{1,n}{1,g},rcv);
                    if flag == 0
                        image{1,n}{1,g} = rcv;
                        continue
                    end
                    % (2)test if the collision point is in the polygon of walls.
                    flag = BoundaryJudge(wall,vertex,point,g);
                    if flag == 0;
                        image{1,n}{1,g} = rcv;
                    end
                end
            end
            
        end
    end
    if n == 2
        for g = 1:1:wnum
            for h = 1:1:wnum
                if h ~= g 
                    image{1,n}{g,h} = image{1,n-1}{1,g} * matrix(:,:,h);
                else
                    image{1,n}{g,h} = rcv;
                end
                
                %calculate the cross point between the trajectory and the
                %relfecting walls
                point = CrossPoint(image{1,n}{g,h},rcv,plane,h);
                % test if the collision point is in the polygon of walls
                flag = BoundaryJudge(wall,vertex,point,h);
                if flag == 1;
                    image{1,n}{g,h} = rcv;
                end
            
                temp = image{1,n}{g,h} * matrix(:,:,h);
                trcv = rcv * matrix(:,:,h);
                %calculate the second order cross point between the
                %trajectory and the relfecting walls
                point2 = CrossPoint(temp,trcv,plane,g);                
                % test if the collision point is in the polygon of walls
                flag = BoundaryJudge(wall,vertex,point2,g);
                if flag == 1;
                    image{1,n}{g,h} = rcv;
                end
                
            end
        end
    end
%     if n == 3
%         for g = 1:1:wnum
%             for k = 1:1:wnum
%                 for h = 1:1:wnum
%                     image{1,n}{(g-1)*wnum+k,h} = image{1,n-1}{g,k} * matrix(:,:,h);
%                 end
%             end
%         end
%     end
end

% compute the IR
distance = cell(1,N);
amplitude = cell(1,N);
time = cell(1,N);
for n = 1:1:wnum
    distance{1,1}{1,n} = norm(rcv-image{1,1}{1,n});
    amplitude{1,1}{1,n} = beta/(4*pi*distance{1,1}{1,n});
    time{1,1}{1,n} = distance{1,1}{1,n}/c;
end
for g = 1:1:wnum
    for h = 1:1:wnum
        distance{1,2}{g,h} = norm(rcv-image{1,2}{g,h});
        amplitude{1,2}{g,h} = beta^2/(4*pi*distance{1,2}{g,h});
        time{1,2}{g,h} = distance{1,2}{g,h}/c;
    end
end

Sample = T*Fs;
TimePoints = ((0:Sample-1)/Fs).';
IR = zeros(Sample,1);

for n = 1:1:wnum
    if time{1,1}{1,n} ~= 0
        IR = IR + amplitude{1,1}{1,n} * sinc((TimePoints-time{1,1}{1,n})*Fs);
    end
end
for g = 1:1:wnum
    for h = 1:1:wnum
        if time{1,2}{g,h} ~= 0
            IR = IR + amplitude{1,2}{g,h} * sinc((TimePoints-time{1,2}{g,h})*Fs);
        end
    end
end

% compute the direct sound power
dd = norm(src-rcv);
da = 1/(4*pi*distance{1,1}{1,n});
dt = dd/c
IR = IR + da * sinc((TimePoints-dt)*Fs);

% plot the IR graph and compare it with real IR
subplot(2,1,1)
plot(TimePoints,IR)

[x,fs] = audioread('00x00y.wav');
subplot(2,1,2)
t = 1:9600+354;
plot(t/Fs,x(t))


% image = src * matrix(:,:,3) * matrix(:,:,1);
% k = image - rcv;
% b = rcv;
% %track = [k(2)+k(3),-k(1),-k(1),-(k(2)*b(1)-k(1)*b(2)+k(3)*b(1)-k(1)*b(3))]
% t = -(dot(plane(1,1:3),b(1:3))+plane(1,4))/dot(plane(1,1:3),k(1:3));
% point = [k(1)*t+b(1),k(2)*t+b(2),k(3)*t+b(3),1]
% 
% image = image * matrix(:,:,1);
% rcv = rcv * matrix(:,:,1);
% k = image - rcv;
% b = rcv;
% %track = [k(2)+k(3),-k(1),-k(1),-(k(2)*b(1)-k(1)*b(2)+k(3)*b(1)-k(1)*b(3))]
% t = -(dot(plane(3,1:3),b(1:3))+plane(3,4))/dot(plane(3,1:3),k(1:3));
% point2 = [k(1)*t+b(1),k(2)*t+b(2),k(3)*t+b(3),1]