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
        cells(k).mesh_size_lower_LG = cell_table.mesh_size_lower_LG{k};
        cells(k).fem_order_lower_LG = cell_table.fem_order_lower_LG{k};
    end

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
        end
        
        % Display status: Cell ID | Progress [Current/Total] | Estimated Time Remaining
        fprintf('Processing cell i=%d | Progress: [%d/%d] | ETR: %s ...\n', ...
            cells(k).i, k, num_cells, etr_str);
        
        % Execute validation for the specific cell
        % (Assumes validate_region_cell is defined in a separate file or below)
        cells(k)
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
    
    total_time = toc(t_start);
    fprintf('Verification Step 2 Completed. Total time: %.2f seconds.\n', total_time);
end