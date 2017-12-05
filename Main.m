%% Initialization
clc; clear

Band_num = 4;
[B,A] = BPFilter(f(Band_num), fs, 10, DispFilter);

Initial_MMR
% Initial_Pollack
% Initial_Tanna
% Initial_Wirth

% choose file type
load (addr1) % type the correct address acording to the Workspace
y = fftdeconv(response, x, 0.99); % compute the measured IR

% choose source and reveiver positions
position = 1; 
% choose the reference direction
direc = 1; 
src = src(position,:); % Speaker postion
rcv = rcv(position,:); % Mic position
ref_direc = ref_direc(direc,:);

% choose display type
overlap_only = 1; % only display overlaped or not 1,0,-1
all = 1; % decide if display the single patha
refl_order = 3; % reflection order 1==first only,2==second only,3==all,-1==none
path = []; % single number represents 1st order reflection, two number vector represents 2nd order reflection
wallrange = -1; % walls that waves will hit in the second reflection -1==none

T = 0.1; % Total time
alpha = 0.001; % air absorbtion coefficient
c = 331.4+0.6*20.7; % velocity of the sound

%% Mirroring processing
N = 2;
image = cell(1,N);
rpoint = cell(1,1);
rpoint2 = cell(1,2);
direction = cell(1,N);
refl_wnum = cell(1,N);
plane = TPlane(wall,wnum,vertex);
for n = 1:1:N
    if n == 1
        for g = 1:1:wnum
            image{1,n}{1,g} = Mirror(wall,vertex,src,g,plane);
            refl_wnum{1,n}{1,g} = g;
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
                refl_wnum{1,n}{1,g} = [];
                rpoint{1,n}{1,g} = [];
            end
            % (2)test if the collision point is in the polygon of walls. if
            % flag = 1, the point is in the polygon. If flag = 0, the point
            % is not in the polygon.
            flag = BoundaryJudge(wall,vertex,point,g);
            if flag == 0;
                image{1,n}{1,g} = rcv;
                refl_wnum{1,n}{1,g} = [];
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
                            refl_wnum{1,n}{1,g} = [];
                            rpoint{1,n}{1,g} = [];
                            disp([g point 1]);
                        end
                    end
                end
                if l ~= g && sum(image{1,n}{1,g}==rcv)~=3 && dot(point-src,plane(l,1:3))~=0
                    tpoint = CrossPoint(point,src,plane,l);
                    temp = BoundaryJudge(wall,vertex,tpoint,l);
                    if temp == 1
                        flag = InPath(tpoint,point,rcv);
                        if flag == 1
                            image{1,n}{1,g} = rcv;
                            refl_wnum{1,n}{1,g} = [];
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
                    refl_wnum{1,n}{g,h} = [g h];
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
                        refl_wnum{1,n}{g,h} = [];
                        rpoint2{1,1}{g,h} = [];
                    end
                    % (2)test if the collision point is in the polygon of walls.
                    flag = BoundaryJudge(wall,vertex,point,g);
                    if flag == 0;
                        image{1,n}{g,h} = rcv;
                        refl_wnum{1,n}{g,h} = [];
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
                                if flag == 1 && ~(abs(sum(tpoint-point))<1e-10)% || abs(sum(tpoint-point2))<1e-10)
                                    image{1,n}{g,h} = rcv;
                                    refl_wnum{1,n}{g,h} = [];
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
                            refl_wnum{1,n}{g,h} = [];
                            rpoint2{1,2}{g,h} = [];
                        end
                        % (1)Test if the point is in the path of the two points
                        flag = InPath(point2,timage,trcv);
                        if flag == 0
                            image{1,n}{g,h} = rcv;
                            refl_wnum{1,n}{g,h} = [];
                            rpoint2{1,2}{g,h} = [];
                        end
                        % (2)test if the collision point is in the polygon of walls
                        flag = BoundaryJudge(wall,vertex,point2,h);
                        if flag == 0;
                            image{1,n}{g,h} = rcv;
                            refl_wnum{1,n}{g,h} = [];
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
                                        refl_wnum{1,n}{g,h} = [];
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
                                        refl_wnum{1,n}{g,h} = [];
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

if refl_order == 1 || refl_order == 3
    for n = 1:1:wnum
        distance{1,1}{1,n} = norm(rcv-image{1,1}{1,n});
        spk_coef = 1;
%         if ~isempty(direction{1,1}{1,n})
%             spk_coef = direc_judge(ref_direc,direction{1,1}{1,n});
%         end
        temp = 1/(4*pi*distance{1,1}{1,n}^2)*exp(-alpha*distance{1,1}{1,n})*spk_coef;
        temp = conv(temp,beta(n,:));
        amplitude{1,1}{1,n} =  temp;
        time{1,1}{1,n} = distance{1,1}{1,n}/c;
    end
end

if refl_order == 2 || refl_order == 3
    for g = 1:1:wnum
        for h = 1:1:wnum
            distance{1,2}{g,h} = norm(rcv-image{1,2}{g,h});
            spk_coef = 1;
            temp = 1/(4*pi*distance{1,2}{g,h}^2)*exp(-alpha*distance{1,2}{g,h})*spk_coef;
            temp = conv(temp,beta(g,:));
            temp = conv(temp,beta(h,:));
            amplitude{1,2}{g,h} =  temp;
            time{1,2}{g,h} = distance{1,2}{g,h}/c;
        end
    end
end

% compute the direct sound power
dd = norm(src-rcv);
da = 1/(4*pi*dd^2);
dt = dd/c;

%% Record sound energy
% Initialization
Sample = ceil((T+dt)*fs);
% TimePoints = ((0:Sample-1)/Fs).';
TimePoints = 0:Sample-1;
IR = zeros(Sample,1);

% Record the direct sound energy
t = TimePoints - round(dt*fs);
t(t==0) = fs*2;
t(t<fs*2) = 0;
t(t==max(t))= 1;
% IR = IR + da * sinc((TimePoints-dt)*Fs*pi);
IR = IR + da * t.';

% Record energy for the first order reflection
if refl_order == 1 || refl_order == 3
    if all == 1;
        for n = 1:1:wnum
            if time{1,1}{1,n} ~= 0 && round(time{1,1}{1,n}*fs)<=Sample
                %           IR = IR + amplitude{1,1}{1,n} * sinc(TimePoints-time{1,1}{1,n}*Fs).';
                t = TimePoints - round(time{1,1}{1,n}*fs);
                t(t==0) = fs*2;
                t(t<fs*2) = 0;
%                 t(t==max(t))= 1;
                t_start = find(t==max(t));
                IR(t_start:t_start+size(amplitude{1,1}{1,n},2)-1) = amplitude{1,1}{1,n};
                disp('Fisrt Reflection');
                disp([n (time{1,1}{1,n}-dt)*1000 max(amplitude{1,1}{1,n})]);
            end
        end
    else
        for n = path:1:path
            if time{1,1}{1,n} ~= 0 && round(time{1,1}{1,n}*fs)<=Sample
                %           IR = IR + amplitude{1,1}{1,n} * sinc(TimePoints-time{1,1}{1,n}*Fs).';
                t = TimePoints - round(time{1,1}{1,n}*fs);
                t(t==0) = fs*2;
                t(t<fs*2) = 0;
%                 t(t==max(t))= 1;
                t_start = find(t==max(t));
                IR(t_start:t_start+size(amplitude{1,1}{1,n},2)-1) = amplitude{1,1}{1,n};
                disp('Fisrt Reflection');
                disp([n (time{1,1}{1,n}-dt)*1000 max(amplitude{1,1}{1,n})]);
            end
        end
    end
end
% Record energy for the second order reflection
if all == 1
    if refl_order == 2 || refl_order == 3
        % wall based ray selection on
        if sum(wallrange)>0
            for uu = 1:1:size(wallrange,2)
                h = wallrange(uu);
                for g = 1:1:wnum
                    if time{1,2}{g,h} ~= 0 && round(time{1,2}{g,h}*fs)<=Sample
                        %               IR = IR + amplitude{1,2}{g,h} * sinc(TimePoints-time{1,2}{g,h}*Fs).';
                        t = TimePoints - round(time{1,2}{g,h}*fs);
                        t(t==0) = fs*2;
                        t(t<fs*2) = 0;
%                         t(t==max(t))= 1;
                        t_start = find(t==max(t));
                        IR(t_start:t_start+size(amplitude{1,2}{g,h},2)-1) = amplitude{1,2}{g,h};
                        disp('Second Reflection');
                        disp([h g (time{1,2}{g,h}-dt)*1000 max(amplitude{1,2}{g,h})]);
                    end
                end
            end
        % wall based ray selection off
        else if wallrange <= 0
                for h = 1:1:wnum
                    for g = 1:1:wnum
                        if time{1,2}{g,h} ~= 0 && round(time{1,2}{g,h}*fs)<=Sample
                            %               IR = IR + amplitude{1,2}{g,h} * sinc(TimePoints-time{1,2}{g,h}*Fs).';
                            t = TimePoints - round(time{1,2}{g,h}*fs);
                            t(t==0) = fs*2;
                            t(t<fs*2) = 0;
%                             t(t==max(t))= 1;
                            t_start = find(t==max(t));
                            IR(t_start:t_start+size(amplitude{1,2}{g,h},2)-1) = amplitude{1,2}{g,h};
                            disp('Second Reflection');
                            disp([h g (time{1,2}{g,h}-dt)*1000 max(amplitude{1,2}{g,h})]);
                        end
                    end
                end
            end
        end
    end
% single path selection on 
else if all == 0
        if time{1,2}{path(2),path(1)} ~= 0 && round(time{1,2}{path(2),path(1)}*fs)<=Sample
            %               IR = IR + amplitude{1,2}{g,h} * sinc(TimePoints-time{1,2}{g,h}*Fs).';
            t = TimePoints - round(time{1,2}{path(2),path(1)}*fs);
            t(t==0) = fs*2;
            t(t<fs*2) = 0;
%             t(t==max(t))= 1;
            t_start = find(t==max(t));
            IR(t_start:t_start+size(amplitude{1,2}{g,h},2)-1) = amplitude{1,2}{g,h};
            disp('Second Reflection');
            disp([path(1) path(2) time{1,2}{path(2),path(1)} max(amplitude{1,2}{path(2),path(1)})]);
        end
    end
end

%% plot the IR graph and compare it with real IR
%     [y,fs] = audioread(diplay_audio);
    [ma,R]=max(y); % Find the direct soudn power, which is the lasrgest one
    IR=IR*ma/max(IR); % normalize the model energy and the measured energy
    R2 = find(IR~=0); % Find the first cell that has value
    di = R-R2(1); % find the position difference between the model and real IR
    
    % Catt acoustic data from Romain
    comp = [9.7 30.17 34.7 35.56 38.36 39.44 42.46 54.53 54.96 59.27 62.28 64.44...
        68.32 71.34 71.98 75 76.72 77.37 78.66 78.88 82.33 87.72 98.06]/1000;
    IR2 = zeros(Sample,1);
    for h = 1:1:size(comp,2)
        t = TimePoints - round(comp(h)*fs);
        t(t==0) = fs*2;
        t(t<fs*2) = 0;
        t(t==max(t))= 1;
        IR2 = IR2 +  0.005*t.';
    end
    IR2(1) = 0.005;
    IR2 = [zeros(R2(1),1);IR2(1:(Sample-R2))];
    
    load('RT.mat');
    IR_t = [IR_t zeros(1,(Sample-size(IR_t,2)))]';
    
    figure;
    if overlap_only == 1
        % match the direct sound position
        if di > 0
            t = 1:Sample-1;
            y = y(di:Sample+di-1);
            hold on
            plot((t)/fs,y(t))
            plot(TimePoints/fs,IR,'LineWidth',2)
%             plot(TimePoints/fs,IR2,'LineWidth',2)
%             plot(TimePoints/fs,IR_t,'LineWidth',2)
%             echo1 = Echo_Density(y',fs);
%             echo2 = Echo_Density(IR,fs);
%             echo3 = Echo_Density(IR2,fs);
%             echo4 = Echo_Density(IR_t,fs);
%             plot(TimePoints/fs, echo1/15)
%             plot(TimePoints/fs, echo2)
%             plot(TimePoints/fs, echo3)
%             plot(TimePoints/fs, echo4)
%             legend('IR','ISM','Catt','RT','Echo-IR','Echo-ISM','Echo-Catt','Echo-RT')
            grid;
            hold off
        else
            t = 1:Sample-1;
            y = [zeros(-di,1);y(1:Sample+di-1)];
            plot((t)/fs,y(t))
            hold on
            plot(TimePoints/fs,IR,'LineWidth',2)
            plot(TimePoints/fs,IR2,'LineWidth',2)
            plot(TimePoints/fs,IR_t,'LineWidth',2)
            echo1 = Echo_Density(y',fs);
            echo1= [zeros(1,round(0.02*fs/2)) echo1];
            echo2 = Echo_Density(IR,fs);
            echo2= [zeros(1,round(0.02*fs/2)) echo2];
            echo3 = Echo_Density(IR2,fs);
            echo3= [zeros(1,round(0.02*fs/2)) echo3];
            echo4 = Echo_Density(IR_t,fs);
            echo4= [zeros(1,round(0.02*fs/2)) echo4];
            plot(0:(Sample-1+round(0.02*fs/2))/fs, echo1/15)
            plot(0:(Sample-1+round(0.02*fs/2))/fs, echo2)
            plot(0:(Sample-1+round(0.02*fs/2))/fs, echo3)
            plot(0:(Sample-1+round(0.02*fs/2))/fs, echo4)
%             plot(TimePoints/fs, echo2)
%             plot(TimePoints/fs, echo3)
%             plot(TimePoints/fs, echo4)
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
                y = y(di:Sample+di);
                plot((t)/fs,y(t))
                hold on
                plot(TimePoints/fs,IR,'LineWidth',2)
                grid;
                hold off
            else
                subplot(3,1,1)
                t = 1:Sample;
                y = [zeros(-di,1);y(1:Sample+di)];
                plot((t)/fs,y(t))
                hold on
                plot(TimePoints/fs,IR,'LineWidth',2)
                grid;
                hold off
            end
            subplot(3,1,2)
            plot((t)/fs,y(t))
            subplot(3,1,3)
            plot(TimePoints/fs,IR,'LineWidth',2)
        end
    end
    
%% visualize the reflection path
% plot the room frame
if refl_order ~= -1
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
    if refl_order == 1 || refl_order == 3
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
        if all == 1
            for u = 1:2:wnum*2-1
                plot3(pathx(:,u:u+1), pathy(:,u:u+1), pathz(:,u:u+1),'Color',color((u+1)/2,:),'LineStyle','-','LineWidth',2)
            end
        else
            for u = path*2-1:2:path*2-1
                plot3(pathx(:,u:u+1), pathy(:,u:u+1), pathz(:,u:u+1),'Color',color((u+1)/2,:),'LineStyle','-','LineWidth',2)
            end
        end
    end
    %plot the second order reflection paths
    if refl_order == 2 || refl_order == 3
        pathx2 = zeros(wnum*3*wnum,2);
        pathy2 = zeros(wnum*3*wnum,2);
        pathz2 = zeros(wnum*3*wnum,2);
        if wallrange ~= -1
            % check if wallrange is correct
            if max(wallrange)>wnum
                error('Wall number required exceeds the max number, please modify the wall range!')
            end
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
        else
            for v = 1:1:wnum
                for u = 1:1:wnum
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
        end
        
        color = zeros(wnum*wnum,3);
        for u = 1:1:wnum*wnum
            % set color arrays. The color gets red when the second
            % value is lower, and it gets blue when the first value is
            % lower, gets yellow when the third value is lower
            %                 color(u,:) = [0.3+0.7*u/(wnum*wnum),0.8*u/(wnum*wnum),0.7+0.2*u/(wnum*wnum)];
            color(u,:) = [0.6+0.4*rand(),0.9*rand(),0.9+0.1*rand()];
        end
        if all == 1
            for u = 1:3:wnum*3*wnum-2
                if sum(pathx2(:,u))~=0 && sum(pathy2(:,u)~=0);
                    plot3(pathx2(:,u:u+2), pathy2(:,u:u+2), pathz2(:,u:u+2),'Color',color((u+2)/3,:),'LineStyle','-','LineWidth',2)
                end
            end
        else
            for u = (path(2)-1)*wnum*3+(path(1)-1)*3+1
                plot3(pathx2(:,u:u+2), pathy2(:,u:u+2), pathz2(:,u:u+2),'Color',color((u+2)/3,:),'LineStyle','-','LineWidth',2)
            end
        end
    end
    grid
    hold off
end
%%