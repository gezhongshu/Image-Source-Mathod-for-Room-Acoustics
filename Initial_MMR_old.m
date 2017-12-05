% MMR structure
dLaser = 0.11;
l = 17+dLaser;
w = 25.57+dLaser;
h = 16.42+dLaser;

% l = 17.1;
% w = 25.5;
% h = 16.5;

% wall = [2 3 7 6; ... % East
%         4 1 5 8; ... % West
%         3 4 8 7; ... % North
%         1 2 6 5; ... % South
%         1 4 3 2; ... % Floor
%         5 6 7 8];    % Ceiling
% 
% vertex = [0 0  0; ...
%         l 0  0; ...
%         l w  0; ...
%         0 w  0; ...
%         0 0  h; ...
%         l 0  h; ...
%         l w  h; ...
%         0 w  h];

wall = [2  3  7  6  0  0  0  0  0  0; ... % East 
        4  1  5  8  0  0  0  0  0  0; ... % West
        3  4  8  7  0  0  0  0  0  0; ... % North
        1  2  6  5  0  0  0  0  0  0; ... % South
        1  4  3  2  0  0  0  0  0  0; ... % Floor
        5  6  7  8  0  0  0  0  0  0; ... % Ceiling
        % The hole at the door
        9  25 26 27 28 29 30 16 24 17; ... % 7
        17 24 23 22 21 20 19 18 0  0; ...  % 8
        26 25 35 31 0  0  0  0  0  0; ...   % 9
        31 32 27 26 0  0  0  0  0  0; ...   % 10
        27 32 33 28 0  0  0  0  0  0; ...   % 11
        33 34 29 28 0  0  0  0  0  0; ...   % 12
        9  17 18 10 0  0  0  0  0  0; ...   % 13
        19 20 12 11 0  0  0  0  0  0; ...   % 14
        20 21 13 12 0  0  0  0  0  0; ...   % 15
        21 22 14 13 0  0  0  0  0  0; ...   % 16
        23 24 16 15 0  0  0  0  0  0; ...   % 17
        % Stairs
        42 43 37 36 0  0  0  0  0  0; ...   % 18
        43 44 38 37 0  0  0  0  0  0; ...   % 19
        44 45 39 38 0  0  0  0  0  0; ...   % 20
        45 46 40 39 0  0  0  0  0  0; ...   % 21
        46 47 41 40 0  0  0  0  0  0; ...   % 22
        ];    

%                             125  250  500  1k   2k   4k   8k
% concrete floor:             0.01 0.03 0.05 0.02 0.02 0.02 0.02
% wood:                       0.09 0.06 0.05 0.05 0.05 0.04 0.04
% Smooth unpainted concrete:  0.01 0.01 0.02 0.02 0.02 0.05 0.05
% Rough lime wash             0.02 0.03 0.04 0.05 0.04 0.03 0.02
concrete = [0.01 0.01 0.02 0.02 0.02 0.05 0.05];
wood = [0.09 0.06 0.05 0.05 0.05 0.04 0.04];
lime = [0.02 0.03 0.04 0.05 0.04 0.03 0.02];
% f = 1; % band number of frequency
abs_coef = [concrete;...     % 1
        concrete;...     % 2
        concrete;...     % 3
        concrete;...     % 4
        wood;...         % 5
        concrete;...     % 6
        % The hole at the door
        lime;...         % 7
        lime;...         % 8
        lime;...         % 9
        concrete;...     % 10
        lime;...         % 11
        lime;...         % 12
        lime;...         % 13
        lime;...         % 14
        lime;...         % 15
        lime;...         % 16
        lime;...         % 17 
        % Stairs
        concrete;...     % 18
        concrete;...     % 19
        concrete;...     % 20
        concrete;...     % 21
        concrete;...     % 22
        ]; % wall absorbtion coefficients
beta = sqrt(1-abs_coef);

vertex = [0 0  0; ...           % 1
        l 0  0; ...             % 2
        l w  0; ...             % 3
        0 w  0; ...             % 4
        0 0  h; ...             % 5
        l 0  h; ...             % 6
        l w  h; ...             % 7
        0 w  h; ...             % 8
        % The hole at the door
        12.08 1.57 0; ...       % 9
        12.08 0    0; ...       % 10
        11.20 0    0; ...       % 11
        8.89  1.45 0; ...       % 12
        5.87  1.45 0; ...       % 13
        3.56  0    0; ...       % 14
        2.67  0    0; ...       % 15
        2.67  1.57 0; ...       % 16
        12.08 1.57 5.13; ...    % 17
        12.08 0    5.13; ...    % 18
        11.20 0    5.13; ...    % 19
        8.89  1.45 5.13; ...    % 20
        5.87  1.45 5.13; ...    % 21
        3.56  0    5.13; ...    % 22 
        2.67  0    5.13; ...    % 23
        2.67  1.57 5.13; ...    % 24
        12.39 1.57 0; ...       % 25
        12.39 1.57 5.76; ...    % 26
        2.36  1.57 5.76; ...    % 27
        2.36  1.57 2.35; ...    % 28
        0     1.57 2.35; ...    % 29
        0     1.57 0; ...       % 30
        12.39 0    5.76; ...    % 31
        2.36  0    5.76; ...    % 32
        2.36  0    2.35; ...    % 33
        0     0    2.35; ...    % 34
        12.39 0    0   ; ...    % 35
        % Stairs
        l     4.42 4.06; ...    % 36
        15.66 4.42 4.06; ...    % 37
        14.20 4.42 3.38; ...    % 38
        12.74 4.42 3.38; ...    % 39
        12.74 2.96 3.38; ...    % 40
        12.74 1.50 2.69; ...    % 41
        l     4.42 5.32; ...    % 42
        15.66 4.42 5.32; ...    % 43
        14.20 4.42 4.64; ...    % 44
        12.74 4.42 4.64; ...    % 45
        12.74 2.96 4.64; ...    % 46
        12.74 1.50 3.95; ...    % 47
        ];

wnum = size(wall,1); % calculate the total wall number
src =[7.53+dLaser, 7.94+dLaser, 3.63+0.15+dLaser; ...
      7.53+dLaser, 7.94+dLaser, 3.63+0.15+dLaser; ...
      7.53+dLaser, 7.94+dLaser, 3.63+0.15+dLaser; ...
      7.53+dLaser, 7.94+dLaser, 3.63+0.15+dLaser
     ];
rcv = [
    8.42+dLaser, 15.14+dLaser, 3.55+0.13+dLaser; ...
    11+dLaser, 17.67-dLaser, 3.55+0.13+dLaser; ...
    13.6780-dLaser, 15.86-dLaser, 3.55+0.13+dLaser; ...
    8.71+dLaser, 16.76-dLaser, 3.55+0.13+dLaser
    ];
ref_direc = [
    0   -1   0; ...
    ];

addr1 = '~/Desktop/reverberation/room-measures/mmr/mmr1response.mat';
addr2 = '~/Desktop/reverberation/room-measures/mmr/mmr2response.mat';
addr3 = '~/Desktop/reverberation/room-measures/mmr/mmr3response.mat';
addr4 = '~/Desktop/reverberation/room-measures/mmr/mmr4response.mat';
addr4b = '~/Desktop/reverberation/room-measures/mmr/mmr4bresponse.mat';
addr4c = '~/Desktop/reverberation/room-measures/mmr/mmr4cresponse.mat';
addr4d = '~/Desktop/reverberation/room-measures/mmr/mmr4dresponse.mat';
addr4e = '~/Desktop/reverberation/room-measures/mmr/mmr4eresponse.mat';