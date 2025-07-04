% This script generates a two-column task list for parallel processing.
% Column 1: Job ID (j)
% Column 2: Status (0 = unprocessed, 1 = completed, 2 = in progress)

fprintf('Generating task list: list_j.csv with initial status 0...\n');

% Define the range of job IDs
job_ids = (1:1220)';  % Create as a column vector

% Create a column for the initial status, all set to 0 (unprocessed)
initial_status = zeros(size(job_ids, 1), 1);

% Combine them into a two-column matrix
output_matrix = [job_ids, initial_status];

% Write the two-column matrix to the CSV file in the 'prep' directory
csvwrite('algo2_list_j.csv', output_matrix);

fprintf('Task list generated successfully.\n');