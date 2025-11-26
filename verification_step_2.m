function verification_step_2(j_list, M_steps, N_steps, output_base_path)
% VERIFICATION_STEP_2 Computes rigorous eigenvalue bounds (Algo 2).
% This function iterates over a grid defined in (x, theta) coordinates.
% It uses Domain Monotonicity to verify the gap lambda_2 < lambda_3.
%
% Inputs:
%   j_list: List of theta indices to process (integer array).
%   M_steps: Number of divisions in x direction (radial).
%   N_steps: Number of divisions in theta direction (angular).
%   output_base_path: Base path string for CSV output (e.g., 'results/step2_bounds').

    if nargin < 2, M_steps = 40; end
    if nargin < 3, N_steps = 200; end
    if nargin < 4, output_base_path = 'results/step2_bounds'; end

    % Define the output file name based on j_list range to support batching
    file_name_result = sprintf('%s_%d_%d.csv', output_base_path, min(j_list), max(j_list));
    
    % Ensure directory exists
    [d, ~, ~] = fileparts(file_name_result);
    if ~isempty(d) && ~exist(d, 'dir'), mkdir(d); end

    % Check if a results file already exists and read its content for resumption
    if isfile(file_name_result)
        existing_data = readmatrix(file_name_result, 'NumHeaderLines', 1);
    else
        existing_data = [];
        % Create new file with header
        headers = ["i", "j", "inf_lam2", "sup_lam2", "inf_lam3", "sup_lam3"];
        writematrix(headers, file_name_result);
    end
    
    % --- Main computation loop ---
    for j = j_list
        fprintf('Processing j = %d / %d ...\n', j, N_steps);
        
        for i = 1:M_steps
            % Check if this cell (i,j) is already computed
            if isempty(existing_data) || ~ismember([i, j], existing_data(:, 1:2), 'rows')

                tic
                % --- 1. Calculate Grid Coordinates (Interval Arithmetic) ---
                % We determine the bounding box for the cell C_{i,j}.
                % Point 1: Lower-Left indices (Smallest domain -> Highest eigenvalues)
                xi_small = get_x_coord(i, M_steps);
                theta_small = get_theta_coord(j, N_steps);
                
                % Point 2: Upper-Right indices (Largest domain -> Lowest eigenvalues)
                xi_large = get_x_coord(i + 1, M_steps);
                theta_large = get_theta_coord(j + 1, N_steps);
                
                % Convert to Cartesian (x, y) for triangle vertex P=(x,y)
                % y = x * tan(theta)
                
                % Smallest Triangle coords (P_small)
                x_s = xi_small;
                y_s = xi_small * tan(theta_small);
                
                % Largest Triangle coords (P_large)
                x_l = xi_large;
                y_l = xi_large * tan(theta_large);

                % --- 2. Define Triangles for FEM ---
                % Vertices: (0,0), (1,0), (x,y)
                
                % 'upper_right_matrix': Represents the geometry that gives the UPPER bound of eigenvalues.
                % By domain monotonicity: Smallest Domain => Largest Eigenvalues.
                % So we use (x_s, y_s).
                upper_right_matrix = [I_intval('0'), I_intval('0'); I_intval('1'), I_intval('0'); x_s, y_s];
                
                % 'lower_left_matrix': Represents the geometry that gives the LOWER bound of eigenvalues.
                % By domain monotonicity: Largest Domain => Smallest Eigenvalues.
                % So we use (x_l, y_l).
                lower_left_matrix = [I_intval('0'), I_intval('0'); I_intval('1'), I_intval('0'); x_l, y_l];
        
                % --- 3. Rigorous Computation ---
                % Filter invalid domains (basic sanity check)
                if sup(x_l) > 0.5
                    
                    % Heuristic for FEM order based on shape thinness (y coordinate)
                    % Thinner triangles require more elements or different strategies.
                    y_check = sup(y_s); 
                    if y_check < 0.2
                        N_rho = 600; fem_ord = 1; N_LG = 0; isLG = 0;
                    elseif y_check < 0.7
                        N_rho = 64; fem_ord = 1; N_LG = 0; isLG = 0;
                    else
                        N_rho = 32; fem_ord = 4; N_LG = 8; isLG = 1;
                    end

                    % Compute bounds
                    % ulams_: Upper bounds of eigenvalues (from Smallest Domain)
                    % returns vector of interval eigenvalues
                    ulams_ = calc_eigen_bounds_any_order(upper_right_matrix, N_LG, N_rho, fem_ord, isLG);
                    
                    % llams_: Lower bounds of eigenvalues (from Largest Domain)
                    llams_ = calc_eigen_bounds_any_order(lower_left_matrix, N_LG, N_rho, fem_ord, isLG);
                    
                    % Compute Hull: [inf(LowerBound), sup(UpperBound)]
                    % blams_(1) is the interval containing the true range of lambda_2
                    % blams_(2) is the interval containing the true range of lambda_3
                    % (Assuming calc_eigen_bounds returns sorted eigenvalues)
                    blams_ = hull(llams_, ulams_);
    
                    % --- 4. Store Results ---
                    % Format: i, j, inf_lam2, sup_lam2, inf_lam3, sup_lam3
                    new_data = [i, j, inf(blams_(1)), sup(blams_(1)), inf(blams_(2)), sup(blams_(2))];
                    writematrix(new_data, file_name_result, 'WriteMode', 'append');
                    
                    % Check Gap: inf(lambda_3) - sup(lambda_2)
                    gap = inf(blams_(2)) - sup(blams_(1));
                    
                    if gap <= 0
                        fprintf('  [WARN] i=%d, j=%d: Gap <= 0 (%.2e)\n', i, j, gap);
                    end
                end
                
                % Optional: print timing per cell
                % t_cost = toc;
                % fprintf('  i=%d done (%.2fs)\n', i, t_cost);
            end
        end
    end
end

% -------------------------------------------------------------------------
% Helper Functions for Coordinate Generation (Interval Arithmetic)
% -------------------------------------------------------------------------

% --- X Coordinate Calculation ---
% Logic: Exponential clustering towards x=0.5
function x_val = get_x_coord(i, M)
    % Constants defined as intervals
    X_START = I_intval('0.5');
    X_END   = I_intval('1.0');
    gamma_x = I_intval('5.0');
    
    % Normalize index u in [0, 1]
    u = (I_intval(i) - I_intval('1')) / I_intval(M);
    
    % Weight function
    weight_x = (exp(gamma_x * u) - I_intval('1')) / (exp(gamma_x) - I_intval('1'));
    
    % Linear mapping based on weight
    x_val = X_START + (X_END - X_START) * weight_x;
end

% --- Theta Coordinate Calculation ---
% Logic: Hybrid clustering (Exponential at both top and bottom edges)
function theta_val = get_theta_coord(j, N)
    % 1. Define Theta Boundaries based on Omega_down definition
    y_min_const = I_intval('0.5') * tan(I_pi/60);
    THETA_GRID_MIN = atan(y_min_const);
    
    val_inside = sqrt(I_intval('3.0')) - I_intval('2.0e-5');
    THETA_GRID_MAX = atan(val_inside);
    
    theta_range = THETA_GRID_MAX - THETA_GRID_MIN;
    
    % 2. Normalize index v in [0, 1]
    v = (I_intval(j) - I_intval('1')) / I_intval(N);
    
    % 3. Hybrid Weighting
    alpha_top = I_intval('8.0');
    alpha_bot = I_intval('3.5');
    
    if v < 0.5
        % Bottom half: Focus on t~0
        v_loc = v * I_intval('2');
        ratio = (exp(alpha_bot * v_loc) - I_intval('1')) / (exp(alpha_bot) - I_intval('1'));
        w_j = I_intval('0.5') * ratio;
    else
        % Top half: Focus on p0
        v_loc = (I_intval('1') - v) * I_intval('2');
        ratio = (exp(alpha_top * v_loc) - I_intval('1')) / (exp(alpha_top) - I_intval('1'));
        w_j = I_intval('1') - I_intval('0.5') * ratio;
    end
    
    % 4. Final Coordinate
    theta_val = THETA_GRID_MIN + w_j * theta_range;
end