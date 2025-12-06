function verification_step_2(input_file,output_file)
    % Define file paths

    % =========================================================================
    % 0. Initialization: Write Header to results.csv
    % =========================================================================
    % This step ensures the file starts with the correct column labels.
    % It overwrites the file if it already exists.
    headers = ["i", "x_inf", "x_sup", "theta_inf", "theta_sup", "lam2_sup", "lam3_inf", "isOK"];
    writematrix(headers, output_file);
    fprintf('Initialized %s with headers.\n', output_file);

    % =========================================================================
    % 1. Read cell_def.csv
    % =========================================================================
    if ~isfile(input_file)
        error('Input file "%s" not found.', input_file);
    end

    % Setup options to read ALL columns as strings
    opts = detectImportOptions(input_file);
    opts.VariableTypes(:) = {'char'}; % Set all column types to string
    cell_table = readtable(input_file, opts);

    % =========================================================================
    % 2. Convert table data to a struct array 'cells'
    % =========================================================================
    num_cells = height(cell_table);
    cells = struct(); 
    
    for k = 1:num_cells
        % Convert strings back to double for numerical processing using str2double.
        % str2double('-') automatically results in NaN, handling the missing values.
        
        cells(k).i = str2double(cell_table.i{k});
        cells(k).x_inf = cell_table.x_inf{k};
        cells(k).x_sup = cell_table.x_sup{k};
        cells(k).theta_inf = cell_table.theta_inf{k};
        cells(k).theta_sup = cell_table.theta_sup{k};
        cells(k).mesh_size_upper = cell_table.mesh_size_upper{k};
        cells(k).fem_order_upper = cell_table.fem_order_upper{k};
        cells(k).mesh_size_lower_cr = cell_table.mesh_size_lower_cr{k};
        cells(k).isLG = str2double(cell_table.isLG{k});
        
        % Handle potential NaNs for LG parameters
        % Even if the CSV had '-', str2double('-') returns NaN correctly.
        cells(k).mesh_size_lower_LG = cell_table.mesh_size_lower_LG(k);
        cells(k).fem_order_lower_LG = cell_table.fem_order_lower_LG(k);
    end

<<<<<<< Updated upstream
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
                    if y_check < 0.1
                        % Very thin triangles
                        N_rho = 600; fem_ord = 1; N_LG = 0; isLG = 0;
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
=======
    % =========================================================================
    % 3. Main processing loop
    % =========================================================================
    t_start = tic; % Start timer for ETR calculation
    
    for k = 1:num_cells
        % --- Progress and ETR Calculation ---
        elapsed_time = toc(t_start);
        if k == 1
            etr_str = "Calculating...";
        else
            avg_time = elapsed_time / (k - 1);       % Average time per cell so far
            remaining_items = num_cells - k + 1;     % Items left (including current)
            est_remaining = avg_time * remaining_items;
            
            % Convert seconds to HH:MM:SS format
            hrs = floor(est_remaining / 3600);
            mins = floor(mod(est_remaining, 3600) / 60);
            secs = round(mod(est_remaining, 60));
            etr_str = sprintf('%02d:%02d:%02d', hrs, mins, secs);
>>>>>>> Stashed changes
        end
        
        % Display status: Cell ID | Progress [Current/Total] | Estimated Time Remaining
        fprintf('Processing cell i=%d | Progress: [%d/%d] | ETR: %s ...\n', ...
            cells(k).i, k, num_cells, etr_str);
        
        % Execute validation for the specific cell
        % (Assumes validate_region_cell is defined in a separate file or below)
        cell_result = validate_region_cell(cells(k));
        
        % Confirm the cell is validated        
        lam2_sup = I_sup(I_intval(char(cell_result.lam2_sup))); 
        lam3_inf = I_inf(I_intval(char(cell_result.lam3_inf))); 
        if lam2_sup<lam3_inf
            isOK ='OK';
        else
            isOK ='NG';
        end
        
        % Prepare data string for CSV output
        % Using string array ensures precise formatting
        str_data = [ ...
            string(cell_result.i), ...
            string(cell_result.x_inf), ...
            string(cell_result.x_sup), ...
            string(cell_result.theta_inf), ...
            string(cell_result.theta_sup), ...
            string(cell_result.lam2_sup), ...
            string(cell_result.lam3_inf), ...
            string(isOK) ...
        ];
        
        % Append result to file
        writematrix(str_data, output_file, 'WriteMode', 'append');
    end
    
<<<<<<< Updated upstream
    % Global counts:
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
    N_steps = length(theta_nodes) - 1;
    
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
=======
    total_time = toc(t_start);
    fprintf('Verification Step 2 Completed. Total time: %.2f seconds.\n', total_time);
>>>>>>> Stashed changes
end