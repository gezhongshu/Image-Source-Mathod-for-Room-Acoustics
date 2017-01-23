%% Initialization
clear;
Initial;

%% generate the general function of every wall
wnum = size(wall,1); % calculate the total wall number
plane = TPlane(wall,wnum,vertex,src);
%% Generate the transforming matrix for each plane
vec1 = [1 0 0]; % reference vector for angular computing
vec2 = [0 1 0]; % reference vector for angular computing
vec3 = [0 0 1]; % vector initial referrence plane xy
Opoint = [0 0 0];
matrix = zeros(4,4,wnum);
% matrix = matrix(plane,wnum,vec,Opoint);
for n = 1:1:wnum
    matrix(:,:,n) = TMatrix(plane,wnum,vec1,vec2,vec3,Opoint,n);
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
        for g = 1:1:1
            % regenerate the matrix by which the image should be multiplied
            % transform the reference point and vector to new coordinate
            n_vec1 = [vec1 1] * matrix(:,:,g);
            n_vec2 = [vec2 1] * matrix(:,:,g);
            n_vec3 = [vec3 1] * matrix(:,:,g); % vector initial referrence plane xy
            n_vec1(4) = [];
            n_vec2(4) = [];
            n_vec3(4) = [];
            n_Opoint = [Opoint 1] * matrix(:,:,g);
            n_Opoint(4) = [];
            n_vec1 = n_vec1 - n_Opoint;
            n_vec2 = n_vec2 - n_Opoint;
            n_vec3 = n_vec3 - n_Opoint;
            % transform the vertex to the new coordinate
            n_vertex = zeros(size(vertex,1),size(vertex,2));
            for h = 1:1:size(vertex,1)
                temp = [vertex(h,:),1] * matrix(:,:,g);
                temp(4) = [];
                n_vertex(h,:) = temp;
            end
            % get the new plane function
            n_plane = TPlane(wall,wnum,n_vertex,image{1,n-1}{1,g});
            % calculate the new matrix
            n_matrix = zeros(4,4,wnum);
            for h = 3%1:1:wnum
                n_matrix(:,:,h) = TMatrix(n_plane,wnum,n_vec1,n_vec2,n_vec3,n_Opoint,h);
            end
            
            for h = 1:1:wnum
                if h ~= g 
                    image{1,n}{g,h} = image{1,n-1}{1,g} * n_matrix(:,:,h);
                else
                    image{1,n}{g,h} = rcv;
                end
                
%                 % Testing of the validity of the point begins
%                 for l = 1:1:wnum
%                     % calculate the cross point between the trajectory and the
%                     % relfecting walls
%                     point = CrossPoint(image{1,n}{g,h},rcv,plane,h);
%                     if l ~= h
%                         % test if there exist obstructions in the path
%                         % (1)Test if the point is in the path of the two points
%                         flag = InPath(point,image{1,n}{g,h},rcv);
%                         if flag == 0
%                             continue
%                         else
%                             % (2)Test if the collision point is in the polygon
%                             % of walls. If the point IS in the polygon, there
%                             % exist obstructions in the path and the point
%                             % should be removed
%                             flag = BoundaryJudge(wall,vertex,point,l);
%                             if flag == 1 && l ~= g
%                                 image{1,n}{g,h} = rcv;
%                             end
%                         end
%                     end
%                     if l == h
%                         % (1)Test if the point is in the path of the two points
%                         flag = InPath(point,image{1,n}{g,h},rcv);
%                         if flag == 0
%                             image{1,n}{g,h} = rcv;
%                             continue
%                         end
%                         % (2)test if the collision point is in the polygon of walls.
%                         flag = BoundaryJudge(wall,vertex,point,h);
%                         if flag == 1;
%                             image{1,n}{g,h} = rcv;
%                         end
%                         temp = image{1,n}{g,h} * matrix(:,:,h);
%                         trcv = rcv * matrix(:,:,h);
%                         %calculate the second order cross point between the
%                         %trajectory and the relfecting walls
%                         point2 = CrossPoint(temp,trcv,plane,g);
%                         % test if the collision point is in the polygon of walls
%                         flag = BoundaryJudge(wall,vertex,point2,g);
%                         if flag == 1;
%                             image{1,n}{g,h} = rcv;
%                         end
%                     end
%                 end

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

%% compute the IR

cTs = c/Fs;

distance = cell(1,N);
amplitude = cell(1,N);
time = cell(1,N);
for n = 1:1:wnum
    distance{1,1}{1,n} = norm(rcv-image{1,1}{1,n});
    amplitude{1,1}{1,n} = beta/(4*pi*distance{1,1}{1,n});
    time{1,1}{1,n} = distance{1,1}{1,n}/c;
end
% for g = 1:1:wnum
%     for h = 1:1:wnum
%         distance{1,2}{g,h} = norm(rcv-image{1,2}{g,h});
%         amplitude{1,2}{g,h} = beta*amplitude{1,1}{1,g};
%         time{1,2}{g,h} = distance{1,2}{g,h}/c;
%     end
% end

Sample = T*Fs;
TimePoints = ((0:Sample-1)/Fs).';
IR = zeros(Sample,1);

for n = 1:1:wnum
    if time{1,1}{1,n} ~= 0
%         IR = IR + amplitude{1,1}{1,n}/cTs * sinc((TimePoints-time{1,1}{1,n})*Fs);
          IR = IR + sinc((TimePoints-time{1,1}{1,n})*Fs);
    end
end
% for g = 1:1:wnum
%     for h = 1:1:wnum
%         if time{1,2}{g,h} ~= 0
% %             IR = IR + amplitude{1,2}{g,h}/cTs * sinc((TimePoints-time{1,2}{g,h})*Fs);
%             IR = IR + sinc((TimePoints-time{1,2}{g,h})*Fs);
%         end
%     end
% end

% compute the direct sound power
dd = norm(src-rcv);
da = 1/(4*pi*dd);
dt = dd/c;
IR = IR +  sinc((TimePoints-dt)*Fs);

% Tw = 2 * round(0.005*Fs);
% cTs = c/Fs;
% Fc = 1;
% 
% distance = cell(1,N);
% amplitude = cell(1,N);
% time = cell(1,N);
% for n = 1:1:wnum
%     distance{1,1}{1,n} = norm(rcv-image{1,1}{1,n});
%     amplitude{1,1}{1,n} = beta/(4*pi*distance{1,1}{1,n});
%     time{1,1}{1,n} = distance{1,1}{1,n};
% end
% for g = 1:1:wnum
%     for h = 1:1:wnum
%         distance{1,2}{g,h} = norm(rcv-image{1,2}{g,h});
%         amplitude{1,2}{g,h} = beta^2/(4*pi*distance{1,2}{g,h});
%         time{1,2}{g,h} = distance{1,2}{g,h};
%     end
% end
% 
% Sample = T*Fs;
% TimePoints = ((0:Sample-1)/Fs).';
% %TimePoints = 1:1:Tw;
% IR = zeros(Sample,1);
% 
% for n = 1:1:wnum
%     if distance{1,1}{1,n} ~= 0
%         fdist = floor(distance{1,1}{1,n});
%         startPosition = floor(distance{1,1}{1,n}/c*Fs-(Tw/2)+1);
%         for l = 1:1:Tw
%             %IR(startPosition+l) = IR(startPosition+l) + amplitude{1,1}{1,n} * 0.5 * (1 - cos(2*pi*((l-(distance{1,1}{1,n}-fdist))/Tw))) * Fc * sinc(pi*Fc*(l-(distance{1,1}{1,n}-fdist)-(Tw/2)));
%              IR(startPosition+l) = IR(startPosition+l) + amplitude{1,1}{1,n} * sinc(pi*Fc*(l-(distance{1,1}{1,n}-fdist)-(Tw/2)));
%         end
%     end
% end
% for g = 1:1:wnum
%     for h = 1:1:wnum
%         if distance{1,2}{g,h} ~= 0
%             fdist = floor(distance{1,2}{g,h});
%             startPosition = floor(distance{1,2}{g,h}/c*Fs-(Tw/2)+1);
%             for l = 1:1:Tw
%                 IR(startPosition+l) = IR(startPosition+l) + amplitude{1,2}{g,h} * 0.5 * (1 - cos(2*pi*((l-(distance{1,2}{g,h}-fdist))/Tw))) * Fc * sinc(pi*Fc*(l-(distance{1,2}{g,h}-fdist)-(Tw/2)));
%             end
%         end
%     end
% end


% plot the IR graph and compare it with real IR
subplot(2,1,1)
plot(TimePoints,IR)

[x,fs] = audioread('00x00y.wav');
subplot(2,1,2)
t = 1:9600+354;
plot(t/Fs,x(t))
%%