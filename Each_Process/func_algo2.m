% --- Main function to compute rigorous eigenvalue bounds over subdomains. ---
% This function iterates over a grid of points (i, j) that parameterize 
% the shape of a triangle. For each small rectangular subdomain R_ij, it computes
% guaranteed upper and lower bounds for the 2nd and 3rd Dirichlet eigenvalues.
% This corresponds to the procedure for the Omega_down^(1) region in the paper.

function func_algo2(j_list)
    
    % Set the output format for numerical display.
    % 'format long infsup' is suitable for interval arithmetic,
    % displaying results as [infimum, supremum].
    
    % Define the output file name for storing results.
    % The name is based on the range of the 'j' index being processed.
    pwd
    file_name_result = ['../results/step2_bounds_' num2str(min(j_list)) '_' num2str(max(j_list)) '.csv']
    
    % Check if a results file already exists and read its content.
    % This prevents re-computation for subdomains that have already been processed.
    if isfile(file_name_result)
        existing_data = readmatrix(file_name_result);
    else
        existing_data = [];
    end
    
    % --- Main computation loop ---
    % Iterates through the grid indices 'j' and 'i' to cover the parameter space.
    for j = j_list
        for i = 1:5000
            % Check if the result for the current grid point (i, j) already exists in the file.
            % If it exists, skip to the next point to avoid redundant calculations.
            if isempty(existing_data) || ~ismember([i, j], existing_data(:, 1:2), 'rows')

                % Get the coordinates of the corners of the current rectangular subdomain R_ij.
                % (x, y) and (tx, ty) define the boundaries of the rectangle in the parameter space.
                % These correspond to the nodes defined in the paper.
                x = xij(i, j); y = yij(i, j);
                tx = xij(i + 1, j); ty = yij(i, j + 1);
        
                % Define the triangle shapes at the corners of the subdomain R_ij.
                % 'upper_right' corresponds to the triangle T^p with p=(x_{i+1}, y_j).
                % 'lower_left' corresponds to the triangle T^p with p=(x_i, y_{j+1}).
                % The vertices are (0,0), (1,0), and the third vertex given by the parameters.
                upper_right = [I_intval('0'), I_intval('0'); I_intval('1'), I_intval('0'); tx, y];
                lower_left = [I_intval('0'), I_intval('0'); I_intval('1'), I_intval('0'); x, ty];
        
                % Proceed with computation only if the triangle vertex is within the moduli space,
                % which is bounded by a unit circle (x^2 + y^2 <= 1).
                if not(tx^2 + y^2 > 1)
                    % Set up the affine transformation matrices S for the perturbation estimate.
                    % These matrices map between triangles at different corners of the subdomain,
                    % used to quantify how eigenvalues change across the subdomain.
                    S_u = [1, (x - tx) / ty; 0, 1];
                    S_l = [1, (tx - x) / y; 0, 1];
                    
                    % Calculate the perturbation factors 'pm' (minimum) and 'pM' (maximum).
                    % These factors provide guaranteed bounds on the change in eigenvalues
                    % across the entire rectangular subdomain R_ij, based on Lemma 4.4.
                    SS_u = S_u^(-1) * (S_u^(-1))';
                    SS_l = S_l^(-1) * (S_l^(-1))';
                    pm = I_intval(I_inf(min(I_veig(SS_u, eye(2, 2), 1))));
                    pM = I_intval(I_sup(max(1 ./ I_veig(eye(2, 2), SS_l, 1))));
                    
                    % Set parameters for the eigenvalue calculation.
                    N_LG = 4; % The Lehmann-Goerisch method is set to compute bounds for the first 4 eigenvalues.
                    fem_ord = 5; % A 5th-degree finite element method is used for high accuracy.

                    % Calculate rigorous eigenvalue bounds at the corners of the subdomain.
                    % 'llams' are bounds for the 'upper_right' triangle vertex (x_{i+1}, y_j).
                    % 'ulams' are bounds for the 'lower_left' triangle vertex (x_i, y_{j+1}).
                    % The function 'isLG' disables the Lehmann-Goerisch method in a specific region
                    % where lambda_3 and lambda_4 are too close to be separated easily.
                    llams = calc_eigen_bounds_any_order(upper_right, N_LG, N_rho(i, j), fem_ord, isLG(x,y));
                    ulams = calc_eigen_bounds_any_order(lower_left, N_LG, N_rho(i, j), fem_ord, isLG(x,y));
    
                    % Apply the perturbation factors to the eigenvalue bounds.
                    % This extends the bounds calculated at the corner points to the entire
                    % rectangular subdomain R_ij, as described in equation (34) of the paper.
                    llams_ = llams * pm;
                    ulams_ = ulams * pM;

                    % Compute the interval hull of the lower and upper estimates.
                    % This gives a single interval [inf, sup] guaranteed to contain the
                    % true eigenvalue for any triangle shape within the subdomain R_ij.
                    blams_ = hull(llams_, ulams_);
    
                    % Store the new results: grid indices (i, j) and the guaranteed interval
                    % bounds for the 2nd (blams_(1)) and 3rd (blams_(2)) eigenvalues.
                    new_data = [i, j, inf(blams_(1)), sup(blams_(1)), inf(blams_(2)), sup(blams_(2))]
                    
                    % Append the newly computed bounds to the CSV file for persistent storage.
                    writematrix(new_data, file_name_result, 'WriteMode', 'append');
                else
                    break;
                end
            else
                break;
            end
        end
    end
end

% --- Function to define the x-coordinate of a grid point (i, j). ---
% The grid is non-uniform; the step size for 'x' changes based on the 'j'
% index. This corresponds to the partitioning of the parameter space Omega_down^(1)
% described in the paper.
function x = xij(i, j)
    if (j >= 1 && j <= 91)
        x = I_intval('0.5') + I_intval('1E-6') * (i - 1);
    elseif (j >= 92 && j <= 181)
        x = I_intval('0.5') + I_intval('1E-5') * (i - 1);
    elseif (j >= 182 && j <= 271)
        x = I_intval('0.5') + I_intval('1E-4') * (i - 1);
    elseif (j >= 272 && j <= 1020)
        x = I_intval('0.5') + I_intval('1E-3') * (i - 1);
    else
        error('Invalid value of j');
    end
end

% --- Function to define the y-coordinate of a grid point (i, j). ---
% The grid is non-uniform; the step size for 'y' changes based on the 'j'
% index. This corresponds to the partitioning of the parameter space Omega_down^(1)
% described in the paper.
function y = yij(i, j)
    sqrt3_over_2 = sqrt(I_intval(3)) / 2;
    if (j >= 1 && j <= 91)
        y = sqrt3_over_2 - I_intval('1E-5') - I_intval('1E-6') * (j - 1);
    elseif (j >= 92 && j <= 181)
        y = sqrt3_over_2 - I_intval('1E-5') - I_intval('1E-6') * 90 - I_intval('1E-5') * (j - 91);
    elseif (j >= 182 && j <= 271)
        y = sqrt3_over_2 - I_intval('1E-5')  - I_intval('1E-6') * 90 - I_intval('1E-5') * 90 - I_intval('1E-4') * (j - 181);
    elseif (j >= 272 && j <= 1020)
        y = sqrt3_over_2 - I_intval('1E-5')  - I_intval('1E-6') * 90 - I_intval('1E-5') * 90 - I_intval('1E-4') * 90 - I_intval('1E-3') * (j - 271);
    else
        error('Invalid value of j');
    end
end

% --- Function to determine the shift parameter 'rho' for the Lehmann-Goerisch method. ---
% 'rho' must be a guaranteed lower bound of an eigenvalue higher than those
% being computed (e.g., rho <= lambda_4 when computing lambda_1 to lambda_3).
% The value is chosen based on the region (defined by 'j').
function val=N_rho(i,j)
    if j<600
        val=32;
    else
        val=50;
    end
end

% --- Function to check if Lehmann-Goerisch (LG) method should be used for lambda_3. ---
% It returns 0 (false) if the point (x,y) is within a small region where 
% lambda_3 and lambda_4 are known to be very close. In this case, a simpler
% bounding method (Lemma 4.2) is used for lambda_3 instead of the LG method
% to avoid difficulties in choosing the shift 'rho'.
% Otherwise, it returns 1 (true).
function val=isLG(x,y)
    if 0.3<y
        val=1;
    else
        val=0;
    end
end