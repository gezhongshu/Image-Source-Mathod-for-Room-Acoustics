% Save the planes structure
clc; clear

% Reference rectangular data
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
src = [4.5, 0.6, 1.9];

% Receiver position
rcv = [1.5, 0.8, 1.9];

