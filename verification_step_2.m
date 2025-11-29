function verification_step_2(j_list, M_steps, output_base_path)
format long infsup

% VERIFICATION_STEP_2 Computes rigorous eigenvalue bounds (Algo 2).
%
% Grid construction logic:
%   - x direction (radial):
%       Determined automatically from {S_i}. A single global M_steps is chosen
%       so that the radial cell size is <= min_i S_i (or some function of S_i).
%   - theta direction (angular):
%       1) Choose coarse breakpoints {theta_i^*} and target diameters {S_i}
%       2) For each coarse interval [theta_i^*, theta_{i+1}^*], define
%          N_i fine angular subdivisions; these N_i are computed from S_i.
%       3) The total number of angular steps is N_steps = sum_i N_i.
%
% Inputs:
%   j_list: List of (global) theta-column indices to process (integer array).
%   output_base_path: Base path string for CSV output.

    if nargin < 2
        output_base_path = 'results/step2_bounds';
    end
    
    % -------------------------------------------------------------
    % Build x and theta grids from coarse angular data {theta_i^*},{S_i}
    % -------------------------------------------------------------
    [x_nodes, theta_nodes, N_steps] = build_global_grids_from_coarse_data(M_steps);
    
    fprintf('Total number of x intervals (N_steps): %d\n', M_steps);
    fprintf('Total number of x nodes: %d\n', M_steps + 1);
    fprintf('Total number of theta intervals (N_steps): %d\n', N_steps);
    fprintf('Total number of theta nodes: %d\n', N_steps + 1);
    
    % If j_list is empty, default to all columns
    if isempty(I_mid(j_list))
        j_list = 1:N_steps;
    end
    
    % Define the output file name based on j_list range to support batching
    file_name_result = sprintf('%s_%d_%d.csv', ...
        output_base_path, I_mid(min(j_list)), I_mid(max(j_list)));
    
    % Ensure directory exists
    [d, ~, ~] = fileparts(file_name_result);
    if ~isempty(d) && ~exist(d, 'dir'), mkdir(d); end

        
    % Check if a results file already exists and read its content for resumption
    if isfile(file_name_result)
        existing_data = readmatrix(file_name_result, 'NumHeaderLines', 1);
        fprintf('Resuming from existing file: %s (%d rows found)\n', ...
                file_name_result, size(existing_data, 1));
    else
        existing_data = [];
        % Create new file with header
        headers = ["i", "j", "inf_lam2", "sup_lam2", "inf_lam3", "sup_lam3"];
        writematrix(headers, file_name_result);
        fprintf('Created new output file: %s\n', file_name_result);
    end
    
    % --- Main computation loop ---
    for j = j_list
        fprintf('========================================\n');
        fprintf('Processing Column j = %d / %d\n', j, N_steps);
        fprintf('========================================\n');
        
        % for i = 1:I_mid(M_steps)
        for i = 1:1
            % Check if this cell (i,j) is already computed
            if isempty(existing_data) || ~ismember([i, j], existing_data(:, 1:2), 'rows')

                tic % Start timer for this cell

                % Linear cell index (for convenience / debugging)
                cell_idx = (j - 1) * M_steps + i;

                % --- 1. Compute probing points for this cell (Interval Arithmetic) ---
                % From cell index -> inner point P_s and outer point P_l.
                [x_s, y_s, x_l, y_l] = get_cell_vertices_from_index( ...
                    cell_idx, x_nodes, theta_nodes);

                % --- [CRITICAL CHECK] Domain Validity ---
                % If the inner point is outside the unit disk, skip this cell.
                dist_sq_s = x_s^2 + y_s^2;
                if I_inf(dist_sq_s) > 1.0
                    continue;
                end
                
                % --- 2. Define Triangles for FEM ---
                upper_right_matrix = [I_intval('0'), I_intval('0'); ...
                                      I_intval('1'), I_intval('0'); ...
                                      x_l,           y_l];
                lower_left_matrix  = [I_intval('0'), I_intval('0'); ...
                                      I_intval('1'), I_intval('0'); ...
                                      x_s,           y_s];

                % --- 3. Rigorous Computation ---
                if I_sup(x_l) > 0.5
                    
                    % Heuristic for FEM order based on triangle thinness
                    y_check = I_sup(y_s)
                    if y_check < 0.2
                        % Very thin triangles
                        N_rho = 1200; fem_ord = 1; N_LG = 0; isLG = 0;
                    elseif y_check < 0.8
                        % Normal triangles
                        N_rho = 32; fem_ord = 1; N_LG = 0; isLG = 0;
                    else
                        % Fat/Equilateral-like triangles
                        % High precision required
                        N_rho = 32; fem_ord = 5; N_LG = 8; isLG = 1;
                    end
                    
                    % Compute bounds
                    ulams_ = calc_eigen_bounds_any_order(lower_left_matrix,  N_LG, N_rho, fem_ord, isLG);
                    llams_ = calc_eigen_bounds_any_order(upper_right_matrix, N_LG, N_rho, fem_ord, isLG);
                    
                    % --- 4. Store Results ---
                    new_data = [i, j, ...
                                I_inf(llams_(1)), I_sup(ulams_(1)), ...
                                I_inf(llams_(2)), I_sup(ulams_(2))];
                    writematrix(new_data, file_name_result, 'WriteMode', 'append');
                    
                    % Check Gap
                    gap = llams_(2) - ulams_(1)
                    t_cost = toc;
                    
                    if gap <= 0
                        fprintf('\n    !!!! WARNING: Gap <= 0 (%.2e) !!!!\n', I_inf(gap));
                        fprintf('    L2 sup: %.6f, L3 inf: %.6f\n', ...
                                I_mid(ulams_(1)), I_mid(llams_(2)));
                    end
                end
            end
        end
    end
    fprintf('Batch Completed.\n');
end


% -------------------------------------------------------------------------
% Build global x and theta nodes from coarse {theta_i^*} and {S_i}
% -------------------------------------------------------------------------
function [x_nodes, theta_nodes, N_steps] = build_global_grids_from_coarse_data(M_steps)
    % 1. Compute theta_min and theta_max from Omega_down definition:
    %    theta_min corresponds to the point
    %       ( sqrt(1 - (0.5*tan(pi/60))^2),  0.5*tan(pi/60) )
    %    theta_max = arctan( sqrt(3) - 2e-5 )

    % y_min = 0.5 * tan(pi/60)
    y_min_const = I_intval('0.5') * tan(I_pi / I_intval('60'));
    % x_min = sqrt(1 - y_min^2)
    x_min_const = sqrt(I_intval('1.0') - y_min_const^2);
    % theta_min = atan(y_min/x_min)
    THETA_MIN = atan(y_min_const / x_min_const);
    
    % theta_max = arctan( sqrt(3) - epsilon ), epsilon = 2e-5
    val_inside = sqrt(I_intval('3.0')) - I_intval('2.0e-5');
    THETA_MAX = atan(val_inside);

    % ---------------------------------------------------------------------
    % 2. Specify coarse angular partition {theta_i^*} and target diameters {S_i}
    % ---------------------------------------------------------------------
    % Number of coarse intervals K (user can change this)
    
    % Coarse breakpoints theta_star_d(i), i = 0,...,K  (user can choose non-uniform)
    theta_star_percent = [I_intval('0'),I_intval('10'),I_intval('95'),I_intval('99.8'),I_intval('99.98'),I_intval('100')];
    theta_star_d = theta_star_percent./I_intval('100').*(THETA_MAX-THETA_MIN)+THETA_MIN;

    % Target diameters S_i (must have length K)
    S_list = [I_intval('1e-3'), I_intval('0.01'), I_intval('1e-3'), I_intval('1e-4'), I_intval('1e-6')];
    
    % Number of coarse intervals K is determined automatically
    K = length(I_mid(theta_star_d)) - 1;     
    
    % ---------------------------------------------------------------------
    % 3. Determine per-band (M_i, N_i) and global M_steps, N_steps
    % ---------------------------------------------------------------------
    % For each coarse interval, we choose:
    %   - N_i: number of angular subdivisions
    %   - M_i: (suggested) number of radial subdivisions
    
    Lx = I_intval('1.0') - I_intval('0.5');
    N_i = I_zeros(1, K);
    M_i = I_zeros(1, K);    
    
    for i = 1:K        
        Ltheta = theta_star_d(i+1) - theta_star_d(i);
        S_i_val = S_list(i);
        if S_i_val <= 0
            error('S_i must be positive.');
        end
        N_i(i) = max(1, ceil(I_sup(Ltheta / S_i_val)));
        M_i(i) = max(1, ceil(I_sup(Lx/ S_i_val)));
    end
    
    % Global counts:
    M_steps = max(M_i);          % choose global radial intervals as max_i M_i
    cumN = cumsum(I_mid(N_i));
    N_steps = cumN(end);
    
    j_start = 1;
    fprintf('=== S_list bands and corresponding j ranges ===\n');
    for b = 1:K
        j_end = cumN(b);
        fprintf('Band %d: S = %.3e,  j = %d .. %d\n', ...
            b, I_mid(S_list(b)), I_mid(j_start), I_mid(j_end));
        j_start = j_end + 1;
    end
    fprintf('==============================================\n');
    
    % ---------------------------------------------------------------------
    % 4. Build theta_nodes (N_steps+1) via per-band uniform subdivision
    % ---------------------------------------------------------------------
    
    % Pre-allocate theta_nodes directly as a column interval vector
    theta_nodes = I_zeros(N_steps + 1, 1);
    
    % Set the first node
    theta_nodes(1) = theta_star_d(1);
    
    current_idx = 1;
    
    for i = 1:K
        a = theta_star_d(i);
        b = theta_star_d(i+1);
        n_loc = N_i(i);          % Interval scalar
        n_val = I_mid(n_loc);    % Integer value for indexing
        
        % Vectorize inner loop: Create m vector [1, 2, ..., n_loc]
        % Note: I_intval(1:n_val) creates a row vector of intervals
        m_vec = I_intval(1:n_val)'; 
        
        % Compute t vector: t = m / n_loc
        t_vec = m_vec / n_loc;
        
        % Compute nodes block: (1 - t)*a + t*b
        % INTLAB supports vector operations for this formula
        nodes_block = (I_intval('1.0') - t_vec) * a + t_vec * b;
        
        % Assign block to global array
        theta_nodes(current_idx + 1 : current_idx + n_val) = nodes_block;
        
        current_idx = current_idx + n_val;
    end
    
    
    % ---------------------------------------------------------------------
    % 5. Build x_nodes (M_steps+1) as uniform on [0.5, 1.0]
    % ---------------------------------------------------------------------
    X_START = I_intval('0.5');
    X_END   = I_intval('1.0');
    
    m_val = I_mid(M_steps);
    
    % Vectorize loop: k = 0:M_steps
    k_vec = I_intval(0:m_val)';  % Column vector
    
    % u = k / M_steps
    u_vec = k_vec / I_intval(M_steps);
    
    % Calculate x_nodes using vector arithmetic
    x_nodes = X_START + (X_END - X_START) * u_vec;
    
end


% -------------------------------------------------------------------------
% Helper: from linear cell index -> all vertices (inner and outer points)
% -------------------------------------------------------------------------
function [x_s, y_s, x_l, y_l] = get_cell_vertices_from_index( ...
    cell_idx, x_nodes, theta_nodes)

    % Infer M_steps and N_steps from node arrays
    M_steps = length(x_nodes)     - 1;
    N_steps = length(theta_nodes) - 1; % Unused but kept for structure

    % Map linear cell index back to (i,j):
    %   cell_idx = (j-1)*M_steps + i
    j = floor((cell_idx - 1) / M_steps) + 1;
    i = cell_idx - (j - 1) * M_steps;
    
    % Column-dependent x subdivision:
    % For each angular column j, subdivide the interval
    % [1/2, sqrt(1 - y_s^2)] into M_steps equal parts, where y_s
    % corresponds to the lower angular boundary theta_nodes(j).
    X_START     = I_intval('0.5');
    theta_small = theta_nodes(j);
    y_for_xmax  = X_START * tan(theta_small);
    x_max       = sqrt(I_intval('1.0') - y_for_xmax^2);
    
    % Inner (small) probing point P_s:
    % uses the lower-left corner of this cell in (x,theta)-grid: (i, j)
    x_s = X_START + (x_max - X_START) * (I_intval(i-1) / I_intval(M_steps));
    y_s = x_s * tan(theta_small);
    
    % Outer (large) probing point P_l:
    % uses the upper-right corner in (x,theta)-grid: (i+1, j+1)
    theta_large = theta_nodes(j+1);
    x_l = X_START + (x_max - X_START) * (I_intval(i) / I_intval(M_steps));
    y_l = x_l * tan(theta_large);
end