%% Initialization

% Initial_Reference_Rectangular;
% src = [2 3.5 2]; % reference rectangular enclosure source position
% rcv = [2 1.5 2]; % reference rectangular enclosure reveiver position

% Initial_Queen_Marry_classroom;
% src = [4.5, 0.6, 1.9]; % Queen Mary classroom source position
% rcv = [1.5, 0.8, 1.9]; % Queen Mary classroom receiver position

% Initial_Queen_Marry_Octagon;
% src = [0, -8, 1.9]; % Queen Mary octagon library
% rcv = [-6, -6, 1.9]; % Queen Mary octagon library

% Initial_Rect;
clc; clear
Initial_MMR

% display preference settings
position = 4; % choose source and reveiver positions
diplay_audio = 'mmr4.wav';
overlap_only = 1; % only display overlaped or not 1,0,-1
imprange = 0; % set the second order impulse that we want to see 0==all -1==none
refl_order = 2; % reflection order 1,2,-1
wallrange = [1:wnum]; % walls that waves will hit in the second reflection

T = 0.3; % Total time
Fs = 48000; % Sample Rate
alpha = 0.001; % air absorbtion coefficient
% beta = sqrt(1-0.05); % Absorbtion Coefficient
c = 331.4+0.6*23.6; % velocity of the sound

if position == 1
    % mmr1
    src = [7.53+dLaser, 7.94+dLaser, 3.63+0.15+dLaser]; % MMR Speaker
    rcv = [8.42+dLaser, 15.14+dLaser, 3.55+0.13+dLaser]; % MMR Mic
else if position == 2
        % mmr2
        src = [7.53+dLaser, 7.94+dLaser, 3.63+0.15+dLaser]; % MMR Speaker
        rcv = [11+dLaser, 17.67-dLaser, 3.55+0.13+dLaser]; % MMR Mic
    else if position == 3
            % mmr3
            src = [7.53+dLaser, 7.94+dLaser, 3.63+0.15+dLaser]; % MMR Speaker
            rcv = [13.6780-dLaser, 15.86-dLaser, 3.55+0.13+dLaser]; % MMR Mic
        else if position == 4
                % mmr4
                src = [7.53+dLaser, 7.94+dLaser, 3.63+0.15+dLaser]; % MMR Speaker
                rcv = [8.71+dLaser, 16.76-dLaser, 3.55+0.13+dLaser]; % MMR Mic
            end
        end
    end
end

%% Mirroring processing
% Wall absorbing coefficients
% Air absorbing coefficients
N = 2;
image = cell(1,N);
rpoint = cell(1,1);
rpoint2 = cell(1,2);
direction = cell(1,N);
plane = TPlane(wall,wnum,vertex);
for n = 1:1:N
    if n == 1
        for g = 1:1:wnum
            image{1,n}{1,g} = Mirror(wall,vertex,src,g,plane);
            %========================================================
            % Testing of the validity of the point begins 
            % Calculate the cross point between the trajectory and the
            % relfecting walls
            point = CrossPoint(image{1,n}{1,g},rcv,plane,g);
            direction{1,n}{1,g} = point - src;
            rpoint{1,n}{1,g} = point;
            % (1)Test if the point is in the path of the two points. If the
            % point is in the path, flag =1, else flag = 0
            flag = InPath(point,image{1,n}{1,g},rcv);
            if flag == 0
                image{1,n}{1,g} = rcv;
                rpoint{1,n}{1,g} = [];
            end
            % (2)test if the collision point is in the polygon of walls. if
            % flag = 1, the point is in the polygon. If flag = 0, the point
            % is not in the polygon.
            flag = BoundaryJudge(wall,vertex,point,g);
            if flag == 0;
                image{1,n}{1,g} = rcv;
                rpoint{1,n}{1,g} = [];
            end
            % (3)test if there exist any obstructions in the path
            for l = 1:1:wnum
                if l ~= g && sum(image{1,n}{1,g}==rcv)~=3 && dot(point-rcv,plane(l,1:3))~=0
                    tpoint = CrossPoint(point,rcv,plane,l);
                    temp = BoundaryJudge(wall,vertex,tpoint,l);
                    if temp == 1
                        flag = InPath(tpoint,point,rcv);
                        if flag == 1
                            image{1,n}{1,g} = rcv;
                            rpoint{1,n}{1,g} = [];
                            disp([g point 1]);
                        end
                    end
                end
            end
            %========================================================
        end
    end
    if n == 2
        for g = 1:1:wnum
            % transform the vertex to the new coordinate
            n_vertex = zeros(size(vertex));
            for h = 1:1:size(vertex,1)
                n_vertex(h,:) = Mirror(wall,vertex,vertex(h,:),g,plane);
            end
            % begin mirroring second order images
            for h = 1:1:wnum
                if h ~= g && sum(image{1,n-1}{1,g}==rcv)~=3
                    n_plane = TPlane(wall,wnum,n_vertex);
                    image{1,n}{g,h} = Mirror(wall,n_vertex,image{1,n-1}{1,g},h,n_plane);
                    %========================================================
                    % Testing of the validity of the point begins
                    % calculate the cross point between the trajectory and the
                    % relfecting walls
                    point = CrossPoint(image{1,n}{g,h},rcv,plane,g);
                    rpoint2{1,1}{g,h} = point;
                    % (1)Test if the point is in the path of the two points
                    flag = InPath(point,image{1,n}{g,h},rcv);
                    if flag == 0
                        image{1,n}{g,h} = rcv;
                        rpoint2{1,1}{g,h} = [];
                    end
                    % (2)test if the collision point is in the polygon of walls.
                    flag = BoundaryJudge(wall,vertex,point,g);
                    if flag == 0;
                        image{1,n}{g,h} = rcv;
                        rpoint2{1,1}{g,h} = [];
                        %                         disp([g h point]);
                    end
                    % (3)test if there exist any obstructions in the path
                    for l = 1:1:wnum
                        % test the first part of the path(from the receiver to the first crosspoint)
                        if l ~= g && sum(image{1,n}{g,h}==rcv)~=3 && dot(point-rcv,plane(l,1:3))~=0
                            tpoint = CrossPoint(point,rcv,plane,l);
                            temp = BoundaryJudge(wall,vertex,tpoint,l);
                            if temp == 1
                                flag = InPath(tpoint,point,rcv);
                                if flag == 1 && ~(abs(sum(tpoint-point))<1e-10 || abs(sum(tpoint-point2))<1e-10)
                                    image{1,n}{g,h} = rcv;
                                    rpoint2{1,1}{g,h} = [];
                                    disp([g h tpoint point]);
                                end
                            end
                        end
                    end
                    
                    % go to the next order reflections
                    if sum(image{1,n}{g,h}==rcv)~=3
                        timage = Mirror(wall,vertex,image{1,n}{g,h},g,plane);
                        trcv = Mirror(wall,vertex,rcv,g,plane);
                        %calculate the second order cross point between the
                        %trajectory and the previous relfecting walls
                        point2 = CrossPoint(timage,trcv,plane,h);
                        direction{1,n}{g,h} = point2 - src;
                        rpoint2{1,2}{g,h} = point2;
                        % test if the point is lying on the cross line of two
                        % planes. If so, one path must be deleted
                        if point(1)==point2(1) && point(2)==point2(2) && point(3)==point2(3) && g > h
                            image{1,n}{g,h} = rcv;
                            rpoint2{1,2}{g,h} = [];
                        end
                        % (1)Test if the point is in the path of the two points
                        flag = InPath(point2,timage,trcv);
                        if flag == 0
                            image{1,n}{g,h} = rcv;
                            rpoint2{1,2}{g,h} = [];
                        end
                        % (2)test if the collision point is in the polygon of walls
                        flag = BoundaryJudge(wall,vertex,point2,h);
                        if flag == 0;
                            image{1,n}{g,h} = rcv;
                            rpoint2{1,2}{g,h} = [];
                            %                         disp([g h point2]);
                        end
                        % (3)test if there exist any obstructions in the path
                        for l = 1:1:wnum
                            % test second part of the path (between two
                            % crosspoints)
                            if l ~= h && sum(image{1,n}{g,h}==rcv)~=3 && dot(point2-point,plane(l,1:3))~=0
                                tpoint = CrossPoint(point2,point,plane,l);
                                temp = BoundaryJudge(wall,vertex,tpoint,l);
                                if temp == 1 
                                    flag = InPath(tpoint,point2,point);
                                    if flag == 1 && ~(abs(sum(tpoint-point))<1e-10 || abs(sum(tpoint-point2))<1e-10)
                                        image{1,n}{g,h} = rcv;
                                        rpoint2{1,2}{g,h} = [];
                                        disp([g h l tpoint]);
                                    end
                                end
                            end
                            % test third part of the path (from source to
                            % the second crosspoint)
                            if l ~= h && sum(image{1,n}{g,h}==rcv)~=3 && dot(point2-src,plane(l,1:3))~=0
                                tpoint = CrossPoint(point2,src,plane,l);
                                temp = BoundaryJudge(wall,vertex,tpoint,l);
                                if temp == 1
                                    flag = InPath(tpoint,point2,src);
                                    if flag == 1 && ~(abs(sum(tpoint-point))<1e-10 || abs(sum(tpoint-point2))<1e-10)
                                        image{1,n}{g,h} = rcv;
                                        rpoint2{1,2}{g,h} = [];
                                        disp([g h tpoint point2]);
                                    end
                                end
                            end
                        end
                    end
                    %========================================================
                else
                    image{1,n}{g,h} = rcv;
                    rpoint2{1,1}{g,h} = [];
                    rpoint2{1,2}{g,h} = [];
                    direction{1,n}{g,h} =[];
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

%% compute the IR
% compute the distances and time spent according to each source position
distance = cell(1,N);
amplitude = cell(1,N);
time = cell(1,N);
for n = 1:1:wnum
    distance{1,1}{1,n} = norm(rcv-image{1,1}{1,n});
    if ~isempty(direction{1,1}{1,n})
        spk_coef = direc_judge(src,rcv,direction{1,1}{1,n});
    end
    amplitude{1,1}{1,n} = beta(n)/(4*pi*distance{1,1}{1,n})*exp(-alpha*distance{1,1}{1,n})*spk_coef;
    time{1,1}{1,n} = distance{1,1}{1,n}/c;
end
for g = 1:1:wnum
    for h = 1:1:wnum
        distance{1,2}{g,h} = norm(rcv-image{1,2}{g,h});
        if ~isempty(direction{1,2}{g,h})
            spk_coef = direc_judge(src,rcv,direction{1,2}{g,h});
        end
        amplitude{1,2}{g,h} = beta(g)*beta(h)/(4*pi*distance{1,2}{g,h})*exp(-alpha*distance{1,2}{g,h})*spk_coef;
        time{1,2}{g,h} = distance{1,2}{g,h}/c;
    end
end

Sample = T*Fs;
% TimePoints = ((0:Sample-1)/Fs).';
TimePoints = 0:Sample-1;
IR = zeros(Sample,1);

% record energy for the first order reflection
for n = 1:1:wnum
    if time{1,1}{1,n} ~= 0 && round(time{1,1}{1,n}*Fs)<=Sample
%           IR = IR + amplitude{1,1}{1,n} * sinc(TimePoints-time{1,1}{1,n}*Fs).';
        t = TimePoints - round(time{1,1}{1,n}*Fs);
        t(t==0) = Fs*2;
        t(t<Fs*2) = 0;
        t(t==max(t))= 1;
        IR = IR + amplitude{1,1}{1,n} * t.';
        disp('Fisrt Reflection');
        disp([n time{1,1}{1,n}]);
    end
end
% record energy for the second order reflection
if imprange ~= -1
    if sum(imprange)>0
        for uu = 1:1:size(imprange,2)
            g = imprange(uu);
            for h = 1:1:wnum
                if time{1,2}{g,h} ~= 0 && round(time{1,2}{g,h}*Fs)<=Sample
                    %               IR = IR + amplitude{1,2}{g,h} * sinc(TimePoints-time{1,2}{g,h}*Fs).';
                    t = TimePoints - round(time{1,2}{g,h}*Fs);
                    t(t==0) = Fs*2;
                    t(t<Fs*2) = 0;
                    t(t==max(t))= 1;
                    IR = IR + amplitude{1,2}{g,h} * t.';
                    disp('Second Reflection');
                    disp([g h time{1,2}{g,h}]);
                end
            end
        end
    else if imprange == 0
            for g = 1:1:wnum
                for h = 1:1:wnum
                    if time{1,2}{g,h} ~= 0 && round(time{1,2}{g,h}*Fs)<=Sample
                        %               IR = IR + amplitude{1,2}{g,h} * sinc(TimePoints-time{1,2}{g,h}*Fs).';
                        t = TimePoints - round(time{1,2}{g,h}*Fs);
                        t(t==0) = Fs*2;
                        t(t<Fs*2) = 0;
                        t(t==max(t))= 1;
                        IR = IR + amplitude{1,2}{g,h} * t.';
                    end
                end
            end
        end
    end
end
% compute the direct sound power
dd = norm(src-rcv);
da = 1/(4*pi*dd);
dt = dd/c;
t = TimePoints - round(dt*Fs);
t(t==0) = Fs*2;
t(t<Fs*2) = 0;
t(t==max(t))= 1;
% IR = IR + da * sinc((TimePoints-dt)*Fs*pi);
IR = IR + da * t.';

%% plot the IR graph and compare it with real IR
if overlap_only == 1 || overlap_only == 0
    [x,fs] = audioread(diplay_audio);
    [ma,R]=max(x); % Find the direct soudn power, which is the lasrgest one
    R2 = find(IR~=0); % Find the first cell that has value
    di = R-R2(1); % find the position difference between the model and real IR
    figure;
    if overlap_only == 1
        % match the direct sound position
        if di > 0
            t = 1:Sample;
            x = x(di:Sample+di);
            plot((t)/Fs,x(t))
            hold on
            plot(TimePoints/Fs,IR,'LineWidth',2)
            grid;
            hold off
        else
            t = 1:Sample;
            x = [zeros(-di,1);x(1:Sample+di)];
            plot((t)/Fs,x(t))
            hold on
            plot(TimePoints/Fs,IR,'LineWidth',2)
            grid;
            hold off
        end
    % plot three graphs including overlapped echograph, real IR and model
    % IR
    else if overlap_only == 0
            % match the direct sound position
            if di > 0
                subplot(3,1,1)
                t = 1:Sample;
                x = x(di:Sample+di);
                plot((t)/Fs,x(t))
                hold on
                plot(TimePoints/Fs,IR,'LineWidth',2)
                grid;
                hold off
            else
                subplot(3,1,1)
                t = 1:Sample;
                x = [zeros(-di,1);x(1:Sample+di)];
                plot((t)/Fs,x(t))
                hold on
                plot(TimePoints/Fs,IR,'LineWidth',2)
                grid;
                hold off
            end
            subplot(3,1,2)
            plot((t)/Fs,x(t))
            subplot(3,1,3)
            plot(TimePoints/Fs,IR,'LineWidth',2)
        end
    end
end

%% visualize the reflection path
% plot the room frame
if refl_order == 1 || refl_order == 2
    figure;
    xx = zeros(wnum*size(wall,2),2);
    yy = zeros(wnum*size(wall,2),2);
    zz = zeros(wnum*size(wall,2),2);
    for u = 1:1:wnum % wall number counter
        for v = 1:1:size(wall,2) % vertice counter for every wall
            if (v == size(wall,2) && wall(u,v)~=0)||(wall(u,v) ~= 0 && wall(u,v+1)==0)
                xx((u-1)*size(wall,2)+v,:) = [vertex(wall(u,v),1),vertex(wall(u,1),1)];
                yy((u-1)*size(wall,2)+v,:) = [vertex(wall(u,v),2),vertex(wall(u,1),2)];
                zz((u-1)*size(wall,2)+v,:) = [vertex(wall(u,v),3),vertex(wall(u,1),3)];
                continue
            else if v < size(wall,2) && wall(u,v+1)~=0
                    xx((u-1)*size(wall,2)+v,:) = [vertex(wall(u,v),1),vertex(wall(u,v+1),1)];
                    yy((u-1)*size(wall,2)+v,:) = [vertex(wall(u,v),2),vertex(wall(u,v+1),2)];
                    zz((u-1)*size(wall,2)+v,:) = [vertex(wall(u,v),3),vertex(wall(u,v+1),3)];
                end
            end
        end
    end
    xx = xx';
    yy = yy';
    zz = zz';
    plot3(xx, yy, zz,'k-','LineWidth',2)
    hold on
    % plot the source and the receiver position
    plot3(src(1),src(2),src(3),'.','MarkerSize',50) % Source position
    plot3(rcv(1),rcv(2),rcv(3),'.','MarkerSize',50) % Receiver position
    % plot the first order reflection paths
    if refl_order == 1
        pathx = zeros(wnum*2,2);
        pathy = zeros(wnum*2,2);
        pathz = zeros(wnum*2,2);
        for u = 1:2:wnum*2-1
            if ~isempty(rpoint{1,1}{1,(u+1)/2})
                pathx(u,:) = [src(1),rpoint{1,1}{1,(u+1)/2}(1)];
                pathx(u+1,:) = [rpoint{1,1}{1,(u+1)/2}(1),rcv(1)];
                pathy(u,:) = [src(2),rpoint{1,1}{1,(u+1)/2}(2)];
                pathy(u+1,:) = [rpoint{1,1}{1,(u+1)/2}(2),rcv(2)];
                pathz(u,:) = [src(3),rpoint{1,1}{1,(u+1)/2}(3)];
                pathz(u+1,:) = [rpoint{1,1}{1,(u+1)/2}(3),rcv(3)];
            end
        end
        pathx = pathx';
        pathy = pathy';
        pathz = pathz';
        color = zeros(wnum,3);
        for u = 1:1:wnum
            % set color arrays. The color gets red when the second
            % value is lower, and it gets blue when the first value is
            % lower, gets yellow when the third value is lower
            %             color(u,:) = [0.3+0.7*u/(wnum),0.8*u/(wnum),0.7+0.2*u/(wnum)];
            color(u,:) = [0.4+0.4*rand(),0.3*rand(),0.9+0.1*rand()];
        end
        for u = 1:2:wnum*2-1
            plot3(pathx(:,u:u+1), pathy(:,u:u+1), pathz(:,u:u+1),'Color',color((u+1)/2,:),'LineStyle','-','LineWidth',2)
        end
    %plot the second order reflection paths
    else if refl_order == 2
            % check if wallrange is correct
            if max(wallrange)>wnum
                error('Wall number required exceeds the max number, please modify the wall range')
            end
            pathx2 = zeros(wnum*3*wnum,2);
            pathy2 = zeros(wnum*3*wnum,2);
            pathz2 = zeros(wnum*3*wnum,2);
            %             for u = 2*3-2:3:2*3-2
            for uu = 1:1:size(wallrange,2)
                u = wallrange(uu);
                for v = 1:1:wnum
                    if ~isempty(rpoint2{1,1}{u,v}) && ~isempty(rpoint2{1,2}{u,v})
                        w = (u-1)*wnum*3+(v-1)*3+1;
                        pathx2(w,:) = [src(1),rpoint2{1,2}{u,v}(1)];
                        pathx2(w+1,:) = [rpoint2{1,2}{u,v}(1),rpoint2{1,1}{u,v}(1)];
                        pathx2(w+2,:) = [rpoint2{1,1}{u,v}(1),rcv(1)];
                        pathy2(w,:) = [src(2),rpoint2{1,2}{u,v}(2)];
                        pathy2(w+1,:) = [rpoint2{1,2}{u,v}(2),rpoint2{1,1}{u,v}(2)];
                        pathy2(w+2,:) = [rpoint2{1,1}{u,v}(2),rcv(2)];
                        pathz2(w,:) = [src(3),rpoint2{1,2}{u,v}(3)];
                        pathz2(w+1,:) = [rpoint2{1,2}{u,v}(3),rpoint2{1,1}{u,v}(3)];
                        pathz2(w+2,:) = [rpoint2{1,1}{u,v}(3),rcv(3)];
                    end
                end
            end
            pathx2 = pathx2';
            pathy2 = pathy2';
            pathz2 = pathz2';
            color = zeros(wnum*wnum,3);
            for u = 1:1:wnum*wnum
                % set color arrays. The color gets red when the second
                % value is lower, and it gets blue when the first value is
                % lower, gets yellow when the third value is lower
%                 color(u,:) = [0.3+0.7*u/(wnum*wnum),0.8*u/(wnum*wnum),0.7+0.2*u/(wnum*wnum)];
                color(u,:) = [0.6+0.4*rand(),0.9*rand(),0.9+0.1*rand()];
            end
            for u = 1:3:wnum*3*wnum-2
                plot3(pathx2(:,u:u+2), pathy2(:,u:u+2), pathz2(:,u:u+2),'Color',color((u+2)/3,:),'LineStyle','-','LineWidth',2)
            end
        end
    end
    grid
    hold off
end
%%