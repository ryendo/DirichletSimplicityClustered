function cell_result = validate_region_cell(cell_data)
    
    % Compute Upper Bound
    tic;
    disp('Upper bound computation');
    cell_ub = cell_upper_eig_bound(cell_data);
    toc;
    
    % Compute Lower Bound
    tic;
    disp('Lower bound computation');
    cell_lb = cell_lower_eig_bound(cell_data);
    toc;
    
    % Validation Logic (Check if Lower(3) > Upper(2))
    disp("Region cell validation:");
    
    % Assuming cell_ub and cell_lb are returned as interval/double arrays
    % where index 2 is lambda_2 and index 3 is lambda_3.
    if cell_lb(3) > cell_ub(2)
        disp("OK");
        % % Calculate relative width for debug/info purposes
        % rel_width = (cell_lb(3) - cell_ub(2)) / (cell_ub(3) - cell_ub(2));
        % fprintf('  Gap Relative Width: %.4e\n', rel_width);
    else
        disp("NG");
        fprintf('  Gap Closed! LB(3)=%.17g, UB(2)=%.17g\n', cell_lb(3), cell_ub(2));
    end

    % Construct the result structure
    cell_result.i = cell_data.i;
    cell_result.x_inf = cell_data.x_inf;
    cell_result.x_sup = cell_data.x_sup;
    cell_result.theta_inf = cell_data.theta_inf;
    cell_result.theta_sup = cell_data.theta_sup;
    
    % Format numerical results to high-precision strings
    cell_result.lam2_sup = compose("%.17g", cell_ub(2));
    cell_result.lam3_inf = compose("%.17g", cell_lb(3));
end