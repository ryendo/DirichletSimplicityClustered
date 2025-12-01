function validate_region_cell(region_cell)
    % region_cell list can be created by
    % [regions,regions_by_idx, xlist,tlist]=get_node_list();
    tic;
    cell_ub = cell_upper_eig_bound(region_cell);
    toc;
    tic;
    cell_lb = cell_lower_eig_bound(region_cell);
    toc;
    
    display("Region cell validation:")
    if cell_lb(3) > cell_ub(2)
        display("OK")
        rel_width = (cell_lb(3) - cell_ub(2))/(cell_ub(3) - cell_ub(2));
        rel_width

	%Task 1: Add code to append [region_cell, cell_lb(2), cell_ub(2),cell_lb(3), cell_up(3)] to file "region_OK.txt"


    else
        display("NG")
	%Task 2:Add code to append [region_cell, cell_lb(2), cell_ub(2),cell_lb(3), cell_up(3)] to file "region_NG.txt"

    end

end
