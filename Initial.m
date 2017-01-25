% Save the planes structure
clc; clear

% wall = [5 6 7 8; ...
%         1 4 8 7; ...
%         1 2 3 4; ...
%         3 2 6 5; ...
%         3 5 8 4; ...
%         1 7 6 2];
% vertex = [0  0  0; ...
%           0  40 0; ...
%           20 40 0; ...
%           20 0  0; ...
%           20 40 25; ...
%           0  40 25; ...
%           0  0  30; ...
%           20 0  30];

% data of Queen Marry's classroom IR
% wall = [1 2 8  7  0 0; ... 
%         2 3 9  8  0 0; ... 
%         3 4 10 9 0 0; ... 
%         4 5 11 10 0 0; ... 
%         5 6 12 11 0 0; ... 
%         1 6 12 7 0 0; ...
%         1 2 3  4  5 6; ...
%         7 8 9 10 11 12];
%     
% vertex = [0  0   0; ...
%           9  0   0; ...
%           9  8   0; ...
%           6  8   0; ...
%           6  7.5 0; ...
%           0  7.5 0; ...
%           0  0   3.5; ...
%           9  0   3.5; ...
%           9  8   3.5; ...
%           6  8   3.5; ...
%           6  7.5 3.5; ...
%           0  7.5 3.5; ...
%           ];

wall = [2 3 7 6; ... % East 
        1 4 8 5; ... % West 
        3 4 8 7; ... % North 
        1 2 6 5; ... % South 
        1 2 3 4; ... % Floor 
        5 6 7 8];    % Ceiling 
vertex = [0 0  0; ...
          5 0  0; ...
          5 4  0; ...
          0 4  0; ...
          0 0  6; ...
          5 0  6; ...
          5 4  6; ...
          0 4  6];

% Source position
% src = [4.5, 0.5, 2.0];
src = [2 3.5 2];

% Receiver position
% rcv = [1.5, 2, 2.0];
rcv = [2 1.5 2];

T = 0.1; % Total time
Fs = 96000; % Sample Rate
beta = 0.8; % Absorbtion Coefficient
c = 343.0; % velocity of the sound
%% Uncomment the following for a 3D plot of the above setup:
 %plot3(SetupStruc.src_traj(:,1),SetupStruc.src_traj(:,2),SetupStruc.src_traj(:,3),'ro-','markersize',4); hold on;
 %plot3(SetupStruc.mic_pos(:,1),SetupStruc.mic_pos(:,2),SetupStruc.mic_pos(:,3),'ko','markerfacecolor',ones(1,3)*.6);
 %axis equal; axis([0 SetupStruc.room(1) 0 SetupStruc.room(2) 0 SetupStruc.room(3)]);
 %box on; xlabel('x-axis (m)'); ylabel('y-axis (m)'); zlabel('z-axis (m)');