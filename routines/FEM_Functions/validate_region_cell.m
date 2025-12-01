function validate_region_cell(region_cell, region_idx)
    % region_cell list can be created by [regions,regions_by_idx, xlist,tlist]=get_node_list();

    h = predict_mesh_size(region_cell);

    tic;
    cell_ub = I_sup(cell_upper_eig_bound(region_cell, 2*h));
    toc;
    tic;
    cell_lb = I_inf(cell_lower_eig_bound(region_cell, h));
    toc;

    fprintf("Region cell validation:\n");

    % Extract values to append
    a2 = cell_lb(2);   % lower bound of 2nd
    b2 = cell_ub(2);   % upper bound of 2nd
    a3 = cell_lb(3);   % lower bound of 3rd
    b3 = cell_ub(3);   % upper bound of 3rd

    % Convert region_cell (which is probably a vector) to row format
    region_vec = region_cell;

    ts = datestr(now, 'yyyy-mm-dd HH:MM:SS.FFF');  % timestamp string
    disp(ts)

    % Concatenate all data into one row
    row = [region_idx, region_vec, a2, b2, a3, b3];
    disp(row)

    if cell_lb(3) > cell_ub(2)
        disp("OK")

        rel_width = (cell_lb(3) - cell_ub(2)) / (cell_ub(3) - cell_ub(2));
        disp(rel_width)

        % =======================
        % Append to region_OK.txt
        % =======================
        fid = fopen("./results/region_OK.txt","a");
        fprintf(fid, '%.17g ', row);
        fprintf(fid, '%s\n', ts);         % append timestamp at end of line
        fclose(fid);

    else
        disp("NG")

        % =======================
        % Append to region_NG.txt
        % =======================
        fid = fopen("./results/region_NG.txt","a");
        fprintf(fid, '%.17g ', row);
        fprintf(fid, '%s\n', ts);         % append timestamp at end of line
        fclose(fid);

    end
end
