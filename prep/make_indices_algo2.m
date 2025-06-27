% MATLAB code: Save a list of numbers from 1 to 1220 to a CSV file.
numbers = 1:1220; % Generate numbers from 1 to 1220.

% ▼▼▼ Change: Save location is now 'prep/list_j.csv' ▼▼▼
% This path is relative to the project root where MATLAB is executed from.
csvwrite('prep/list_j.csv', numbers');