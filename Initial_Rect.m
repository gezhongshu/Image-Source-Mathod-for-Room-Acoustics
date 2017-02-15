% Save the planes structure
clc; clear

wall = [2 3 7 6; ... % East
    1 4 8 5; ... % West
    3 4 8 7; ... % North
    1 2 6 5; ... % South
    1 2 3 4; ... % Floor
    5 6 7 8];    % Ceiling

l = 6.5;
w = 4.2;
h = 2.7;

vertex = [0 0  0; ...
    l 0  0; ...
    l 6.5  0; ...
    0 6.5  0; ...
    0 0  2.7; ...
    l 0  2.7; ...
    l 6.5  2.7; ...
    0 6.5  2.7];



% Source position
src = [5.30, 2.10, 0.78];

% Receiver position
rcv = [1.50, 3.10, 1.65];

