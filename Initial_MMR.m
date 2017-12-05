% Load the measured source and its sample rate
load ('~/Desktop/reverberation/room-measures/source.mat')
% Load measured IR
addr1 = '~/Desktop/reverberation/room-measures/mmr/mmr1response.mat';
addr2 = '~/Desktop/reverberation/room-measures/mmr/mmr2response.mat';
addr3 = '~/Desktop/reverberation/room-measures/mmr/mmr3response.mat';
addr4 = '~/Desktop/reverberation/room-measures/mmr/mmr4response.mat';
addr4b = '~/Desktop/reverberation/room-measures/mmr/mmr4bresponse.mat';
addr4c = '~/Desktop/reverberation/room-measures/mmr/mmr4cresponse.mat';
addr4d = '~/Desktop/reverberation/room-measures/mmr/mmr4dresponse.mat';
addr4e = '~/Desktop/reverberation/room-measures/mmr/mmr4eresponse.mat';

% Space volume
V = 7.2630e+03;

% Different source and receiver positions in the model
dLaser = 0.11;
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

% Wall absorption coefficients at different frequency bands
f = [125 250 500 1000 2000 4000 8000];
cinderblock = [0.02 0.02 0.03 0.03 0.04 0.05 0.05];
wood = [0.42 0.21 0.1 0.08 0.06 0.06 0.06];
plaster = [0.02 0.02 0.03 0.03 0.04 0.05 0.05];
steel_frame = [0.15 0.1 0.06 0.04 0.04 0.05 0.05];

% f_gypsum = Wall_Filter(f, gypsum, fs,  DispFilter);
% f_wood = Wall_Filter(f, wood, fs, DispFilter);
% f_plaster = Wall_Filter(f, plaster, fs, DispFilter);
% f_steel_frame = Wall_Filter(f, steel_frame, fs, DispFilter);

f_cinderblock = sqrt(1-cinderblock(Band_num));
f_wood = sqrt(1-wood(Band_num));
f_plaster = sqrt(1-plaster(Band_num));
f_steel_frame = sqrt(1-steel_frame(Band_num));

beta = [f_cinderblock;...     % 1
    f_cinderblock;...     % 2
    f_cinderblock;...     % 3
    f_cinderblock;...     % 4
    f_wood;...         % 5
    f_steel_frame;...     % 6
%     The hole at the door
    f_plaster;...         % 7
    f_plaster;...         % 8
    f_plaster;...         % 9
    f_cinderblock;...     % 10
    f_plaster;...         % 11
    f_cinderblock;...         % 12
    f_plaster;...         % 13
    f_cinderblock;...         % 14
    f_cinderblock;...         % 15
    f_cinderblock;...         % 16
    f_plaster;...         % 17
%     Stairs
    f_cinderblock;...     % 18
    f_cinderblock;...     % 19
    f_cinderblock;...     % 20
    f_cinderblock;...     % 21
    f_cinderblock;...     % 22
    ]; % wall absorbtion coefficients

% MMR structure
l = 17+dLaser;
w = 25.57+dLaser;
h = 16.42+dLaser;

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
%     42 43 37 36 0  0  0  0  0  0; ...   % 18
%     43 44 38 37 0  0  0  0  0  0; ...   % 19
%     44 45 39 38 0  0  0  0  0  0; ...   % 20
%     45 46 40 39 0  0  0  0  0  0; ...   % 21
%     46 47 41 40 0  0  0  0  0  0; ...   % 22
    ];

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