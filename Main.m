%% Initialization
clear;
Initial;

%% Mirroring processing
% Check the vailidity of the path
% Wall absorbing coefficients
% Air absorbing coefficients
wnum = size(wall,1); % calculate the total wall number
N = 2;
image = cell(1,N);
plane = TPlane(wall,wnum,vertex,src);
for n = 1:1:N
    if n == 1
        for g = 1:1:wnum
            image{1,n}{1,g} = Mirror(wall,vertex,src,g);
            %========================================================
            % Testing of the validity of the point begins calculate the
            % cross point between the trajectory and the relfecting
            % walls
            point = CrossPoint(image{1,n}{1,g},rcv,plane,g);
            % (1)Test if the point is in the path of the two points
            flag = InPath(point,image{1,n}{1,g},rcv);
            if flag == 0
                image{1,n}{1,g} = rcv;
            end
            % (2)test if the collision point is in the polygon of walls.
            flag = BoundaryJudge(wall,vertex,point,g);
            if flag == 0; 
                image{1,n}{1,g} = rcv;
            end
            % (3)test if there exist any obstructions in the path
            for l = 1:1:wnum
                if l ~= g && sum(image{1,n}{1,g}==rcv)==0
                    tpoint = CrossPoint(point,rcv,plane,l);
                    flag = InPath(tpoint,point,rcv);
                    if flag == 1
                        image{1,n}{1,g} = rcv;
                        disp([g point]);
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
                n_vertex(h,:) = Mirror(wall,vertex,vertex(h,:),g);
            end
            % begin mirroring second order images
            for h = 1:1:wnum
                if h ~= g %&& sum(image{1,n}{1,g}==rcv)==0
                    image{1,n}{g,h} = Mirror(wall,n_vertex,image{1,n-1}{1,g},h);
                    %========================================================
                    % Testing of the validity of the point begins
                    % calculate the cross point between the trajectory and the
                    % relfecting walls
                    point = CrossPoint(image{1,n}{g,h},rcv,plane,g);
                    % (1)Test if the point is in the path of the two points
                    flag = InPath(point,image{1,n}{g,h},rcv);
                    if flag == 0
                        image{1,n}{g,h} = rcv;
                    end
                    % (2)test if the collision point is in the polygon of walls.
                    flag = BoundaryJudge(wall,vertex,point,g);
                    if flag == 0;
                        image{1,n}{g,h} = rcv;
%                         disp([g h point]);
                    end
                    % (3)test if there exist any obstructions in the path
                    for l = 1:1:wnum
                        if l ~= g && sum(image{1,n}{g,h}==rcv)==0
                            tpoint = CrossPoint(point,rcv,plane,l);
                            flag = InPath(tpoint,point,rcv);
                            if flag == 1
                                image{1,n}{g,h} = rcv;
                                disp([g h tpoint point]);
                            end
                        end
                    end
                    
                    % go to the next order reflections
                    timage = Mirror(wall,vertex,image{1,n}{g,h},g);
                    trcv = Mirror(wall,vertex,rcv,g);
                    %calculate the second order cross point between the
                    %trajectory and the previous relfecting walls
                    point2 = CrossPoint(timage,trcv,plane,h);
                    % test if the point is lying on the cross line of two
                    % planes, if so, one path must be deleted
                    if point(1)==point2(1) && point(2)==point2(2) && point(3)==point2(3) && g > h
                        image{1,n}{g,h} = rcv;
                    end
                    % (1)Test if the point is in the path of the two points
                    flag = InPath(point2,timage,trcv);
                    if flag == 0
                        image{1,n}{g,h} = rcv;
                    end
                    % (2)test if the collision point is in the polygon of walls
                    flag = BoundaryJudge(wall,vertex,point2,h);
                    if flag == 0;
                        image{1,n}{g,h} = rcv;
                    end
                    % (3)test if there exist any obstructions in the path
                    for l = 1:1:wnum
                        if l ~= h && image{1,n}{g,h}(1)~=rcv(1) && image{1,n}{g,h}(2)~=rcv(2) && image{1,n}{g,h}(3)~=rcv(3)
                            tpoint = CrossPoint(point2,rcv,plane,l);
                            flag = InPath(tpoint,point2,rcv);
                            if flag == 1
                                image{1,n}{g,h} = rcv;
                                disp([g h tpoint point2]);
                            end
                        end
                    end
                    %========================================================
                else
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

%% compute the IR

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
        amplitude{1,2}{g,h} = beta*beta/(4*pi*distance{1,2}{g,h});
        time{1,2}{g,h} = distance{1,2}{g,h}/c;
    end
end

Sample = T*Fs;
% TimePoints = ((0:Sample-1)/Fs).';
TimePoints = 0:Sample-1;
IR = zeros(Sample,1);

for n = 1:1:wnum
    if time{1,1}{1,n} ~= 0
%         IR = IR + amplitude{1,1}{1,n} * sinc((TimePoints-time{1,1}{1,n})*Fs*pi);
        t = TimePoints - round(time{1,1}{1,n}*Fs);
        t(t==0) = Fs*2;
        t(t<Fs*2) = 0;
        t(t==max(t))= 1;
        IR = IR + amplitude{1,1}{1,n} * t.';
    end
end
for g = 1:1:wnum
    for h = 1:1:wnum
        if time{1,2}{g,h} ~= 0
%             IR = IR + amplitude{1,2}{g,h} * sinc((TimePoints-time{1,2}{g,h})*Fs*pi);
            t = TimePoints - round(time{1,2}{g,h}*Fs);
            t(t==0) = Fs*2;
            t(t<Fs*2) = 0;
            t(t==max(t))= 1;
            IR = IR + amplitude{1,2}{g,h} * t.';
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
subplot(2,1,1)
plot(TimePoints,IR)

[x,fs] = audioread('00x00y.wav');
subplot(2,1,2)
t = 1:9600+354;
plot(t/Fs,x(t))
%%