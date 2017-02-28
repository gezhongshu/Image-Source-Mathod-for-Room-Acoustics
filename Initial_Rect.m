% Save the planes structure
clc; clear

wall = [2 3 7 6; ... % East
    4 1 5 8; ... % West
    3 4 8 7; ... % North
    1 2 6 5; ... % South
    1 4 3 2; ... % Floor
    5 6 7 8];    % Ceiling

l = 6.5;
w = 4.2;
h = 2.7;

vertex = [0 0  0; ...
        l 0  0; ...
        l w  0; ...
        0 w  0; ...
        0 0  h; ...
        l 0  h; ...
        l w  h; ...
        0 w  h];



% Source position
src = [5.30, 2.10, 0.78];

% Receiver position
rcv = [1.50, 3.10, 1.65];

