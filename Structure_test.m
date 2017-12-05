% Initial_MMR
% Initial_Pollack
% Initial_Tanna
Initial_Wirth

% Range Selection
num1 = 1;
num2 = wnum;

figure;
    xx = zeros(wnum*size(wall,2),2);
    yy = zeros(wnum*size(wall,2),2);
    zz = zeros(wnum*size(wall,2),2);
    for u = num1:1:num2 % wall number counter
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
    
    plane = TPlane(wall,wnum,vertex);
    
    src = src(1,:); % Speaker postion
%     rcv = rcv(1,:); % Mic position
    plot3(src(1),src(2),src(3),'.','MarkerSize',50) % Source position
    plot3(rcv(1,1),rcv(1,2),rcv(1,3),'.','MarkerSize',50) % Receiver position
    plot3(rcv(2,1),rcv(2,2),rcv(2,3),'.','MarkerSize',50) % Receiver position
    plot3(rcv(3,1),rcv(3,2),rcv(3,3),'.','MarkerSize',50) % Receiver position

%     for u = num1:1:num2 % wall number counter
%         plane(u,1:3) = plane(u,1:3)/norm(plane(u,1:3));
%         temp = vertex(wall(u,1),:) + plane(u,1:3);
%         if u == 4
%             disp([plane(u,:)])
%         end
%         plot3([vertex(wall(u,1),1) temp(1)],[vertex(wall(u,1),2) temp(2)],[vertex(wall(u,1),3) temp(3)],'Color',[0 0.5 0],'LineStyle','-','LineWidth',2)
%     end