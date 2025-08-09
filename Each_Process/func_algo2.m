% --- Main function to compute rigorous eigenvalue bounds using domain monotonicity. ---
% This function iterates over a grid defined in (r, h) coordinates, where
% r = sqrt(x^2+y^2) and h = y. It computes guaranteed bounds for eigenvalues
% in the Omega_down region, where eigenvalues are well-separated.
% Instead of perturbation theory, it leverages the domain monotonicity property of
% Dirichlet eigenvalues (Lemma from Fig. 9, Eq. 39 in the paper).

function func_algo2(j_list)
    
    % Define the output file name for storing results for the Omega_down region.
    file_name_result = ['../results/step2_bounds_' num2str(min(j_list)) '_' num2str(max(j_list)) '.csv'];
    
    % Check if a results file already exists and read its content to avoid re-computation.
    if isfile(file_name_result)
        existing_data = readmatrix(file_name_result);
    else
        existing_data = [];
    end
    
    % --- Main computation loop ---
    % Iterates through the grid indices 'j' and 'i' to cover the parameter space.
    for j = j_list
        for i=1:20 %LIU
            tic
            % Check if the result for the current grid cell (i, j) already exists.
            if isempty(existing_data) || ~ismember([i, j], existing_data(:, 1:2), 'rows')

                % Get the (r, h) coordinates for the corners of the current grid cell.
                % The cell is defined by [ri(i), ri(i+1)] x [hj(j), hj(j+1)].
                r = ri(i); h = hj(j);
                r_ = ri(i+1); h_ = hj(j + 1);
                
                % Convert (r, h) coordinates back to Cartesian (x, y) coordinates
                % to define the triangle vertices.
                x = sqrt(r^2-h^2); y = h;
                x_ = sqrt(r_^2-h_^2); y_ = h_;

                disp('Value of r,h')
                [r,h]
                disp('Value of r_,h_')
                [r_,h_]
                disp('Value of x,y')
                [x,y]
                disp('Value of x_,y_')
                [x_,y_]

                % --- Define triangles at the corners of the (r,h) grid cell ---
                % According to domain monotonicity, the smallest triangle (at r_{i+1}, h_{j+1})
                % gives the highest eigenvalues, which serve as the lower bound for the cell.
                % The largest triangle (at r_i, h_j) gives the lowest eigenvalues, which serve
                % as the upper bound for the cell.
                % 'upper_right' corresponds to the smallest triangle in the cell.
                % 'lower_left' corresponds to the largest triangle in the cell.
                upper_right = [I_intval('0'), I_intval('0'); I_intval('1'), I_intval('0'); x_, y_];
                lower_left = [I_intval('0'), I_intval('0'); I_intval('1'), I_intval('0'); x, y];
        
                % Proceed only if the point is within the valid parameter space for Omega_down^(2).
                if r<1 && not(y_<tan(I_pi/60) / 2) && x_>0.5 && r>h
                    % Set parameters for the eigenvalue calculation.
                    % N_rho: A rough shift for the bounding method.
                    % fem_ord = 1: 1st-degree (linear) finite elements are sufficient here as
                    % eigenvalues are well-separated.
                    
                    % Calculate eigenvalue bounds using the basic projection method (Lemma 4.2).
                    % The Lehmann-Goerisch method is disabled (second argument is 0).
                    % 'llams_': Lower eigenvalue bounds for the cell, computed on the smallest triangle.
                    % 'ulams_': Upper eigenvalue bounds for the cell, computed on the largest triangle.

                    if upper_right(6) < 0.2
                        N_rho = 600; fem_ord = 1; N_LG = 0; isLG = 0;
                    elseif upper_right(6) < 0.7
                        N_rho = 64; fem_ord = 1; N_LG = 0; isLG = 0;
                    else
                        N_rho = 32; fem_ord = 4; N_LG = 8; isLG = 1;
                    end

                    llams_ = calc_eigen_bounds_any_order(upper_right, N_LG, N_rho, fem_ord, isLG);
                    ulams_ = calc_eigen_bounds_any_order(lower_left, N_LG, N_rho, fem_ord, isLG);
                    
                    % Compute the interval hull to get the final guaranteed bounds for the entire cell.
                    % This implements the domain monotonicity principle from Eq. (35).
                    blams_ = hull(llams_, ulams_)
    
                    % Store the new results: grid indices and guaranteed bounds.
                    new_data = [i, j, inf(blams_(1)), sup(blams_(1)), inf(blams_(2)), sup(blams_(2))]
                    % Append the newly computed bounds to the CSV file.
                    writematrix(new_data, file_name_result, 'WriteMode', 'append');
                    if sup(blams_(1)) < inf(blams_(2))
                        disp ("Verification Successful")
                    else
                        disp ("Unsuccessful!!!!!!!")
                    end
                    toc
                end
            else
                break;
            end
        end
    end
end

% --- Function to define the r-coordinate of a grid point (i). ---
% r = sqrt(x^2 + y^2). The grid points are defined for the domain
% monotonicity approach.
function ri_val = ri(i)
    ri_val = sqrt(I_intval('0.5')^2 + (tan(I_pi/60)/2)^2) + (1 - sqrt(I_intval('0.5')^2 + (tan(I_pi/60)/2)^2)) * (i-1) / I_intval('19');
end

function hj_val = hj(j)
    if j >= 1 && j <= 100
        % Explicit formula for the first half of the domain (j=1 to 100)
        hj_val = tan(I_pi/60)/2 + (sqrt(I_intval(3))/2 - I_intval('1e-5') - tan(I_pi/60)/2) * I_intval('0.5') * (2*(j-1)/I_intval('199'))^I_intval('1.5');
    else
        % Explicit formula for the second half of the domain (j=101 to 200)
        hj_val = tan(I_pi/60)/2 + (sqrt(I_intval(3))/2 - I_intval('1e-5') - tan(I_pi/60)/2) * (1 - I_intval('0.5') * (2*(1 - (j-1)/I_intval('199')))^I_intval('1.5'));
    end
end