function status = validate_region_cell_auto_adjust(region_cell, cell_idx)
    % Ensure cell_idx is always a char row vector
    cell_idx = char(cell_idx);

    % region_cell list can be created by:
    % [regions,regions_by_idx, xlist,tlist] = get_node_list();

    % cell_idx format example: "123.1.2.2.4"
    parts = strsplit(cell_idx, '.');   % use strsplit for char
    level_number = numel(parts);

    if level_number > 5
        disp("Maximum subdivision level reached. Marking as NG.");
        status = 0;
        return;
    end

    global region_bound_validation_file_OK region_bound_validation_file_NG
    disp("Region Cell Index with auto size adjustment:")
    disp(cell_idx)

    predicted_h = predict_mesh_size(region_cell);
    if predicted_h < 0.002
        disp("predicted mesh size:")
        disp(predicted_h)
    end
    h = max(0.002, predicted_h);

    disp("mesh size:")
    disp(h)

    % Compute bounds
    tic; cell_ub = I_sup(cell_upper_eig_bound(region_cell, h)); toc;
    tic; cell_lb = I_inf(cell_lower_eig_bound(region_cell, h)); toc;

    % ==========================================================
    %        Decide subdivision vs. success (same as before)
    % ==========================================================
    if cell_lb(3) <= cell_ub(2)
        status = 0;

        mid_x = 0.5 * (region_cell(1) + region_cell(2));
        mid_t = 0.5 * (region_cell(3) + region_cell(4));

        sub_cell1 = [region_cell(1), mid_x, region_cell(3), mid_t];
        sub_cell2 = [mid_x, region_cell(2), region_cell(3), mid_t];
        sub_cell3 = [region_cell(1), mid_x, mid_t, region_cell(4)];
        sub_cell4 = [mid_x, region_cell(2), mid_t, region_cell(4)];

        sub_cells = {sub_cell1, sub_cell2, sub_cell3, sub_cell4};
        sub_cell_status = 1;

        idx = 1;
        for sc = sub_cells
            % Construct next index safely as char
            next_idx = [cell_idx, '.', num2str(idx)];
            status_sub = validate_region_cell_auto_adjust(sc{1}, next_idx);

            if status_sub == 0
                sub_cell_status = 0;
            end
            idx = idx + 1;
        end

        if sub_cell_status == 1
            status = 1;
            return  % All sub-cells validated successfully, not need to output log anymore.
        end

    else
        status = 1;
    end


    % ==========================================================
    %     Output logging for OK / NG (same behavior as before)
    % ==========================================================
    fprintf("Region cell validation:\n");

    lb2 = cell_lb(2);
    up2 = cell_ub(2);
    lb3 = cell_lb(3);
    ub3 = cell_ub(3);

    region_vec = region_cell;
    ts = datestr(now, 'yyyy-mm-dd HH:MM:SS.FFF');

    row = [region_vec, lb2, up2, lb3, ub3];
    disp(ts)
    disp(row)

    if status == 1
        disp("OK")

        rel_width = (cell_lb(3) - cell_ub(2)) / (cell_ub(3) - cell_ub(2));
        disp("Relative width of lower(3) - upper(2):")
        disp(rel_width)

        row_2 = [rel_width, h];
        fid = fopen(region_bound_validation_file_OK, "a");
    else
        disp("NG")
        row_2 = [-1, h];
        fid = fopen(region_bound_validation_file_NG, "a");
    end

    % Write log line — cell_idx is now guaranteed a char
    fprintf(fid, '%s ', cell_idx);
    fprintf(fid, '%.17g ', row);
    fprintf(fid, '%.6f ', row_2);
    fprintf(fid, '%s\n', ts);
    fclose(fid);

end