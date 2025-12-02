function [cell_list, xlist, tlist]=get_node_list()

    %from left to right
    cell_list_row = get_cells_for_same_y(0.5); % x1,x2,t1,t2;
    last_cell = cell_list_row(end,:);

    y = last_cell(1)*tan(last_cell(4));

    N = size(cell_list_row,1);
    tlist = zeros(N+1,1);  
    t_idx = 1;
    tlist(t_idx) = last_cell(3);

    xlist = zeros(N+1,1);
    ylist = zeros(N+1,1);

    cell_list =[];
    
    for k = N:-1:1
        current_cell = cell_list_row(k,:);
        t1 = current_cell(3);
        t2 = current_cell(4);
        theta = t2;

        t_idx = t_idx +1;
        tlist(t_idx) = t2;  

        for idx = k:N
            x1 = cell_list_row(idx,1);
            x2 = cell_list_row(idx,2);
            cell_list = [cell_list; x1,x2,t1,t2];
            if k==1
                y1 = tan(t2)*x1;
                y2 = tan(t2)*x2;
                xlist(idx) = x1; ylist(idx) = y1;
            end
        end
    end
    xlist(N+1) = x2; ylist(N+1) = y2;

    y = ylist(1);
    x = xlist(1);

    k=1;
    while y < 0.1
        y
        [new_xlist,new_theta] = get_x_nodex_for_same_theta(xlist,theta);
        xlist = new_xlist;
        local_cells = [xlist(1:end-1),xlist(2:end),theta*ones(N,1),new_theta*ones(N,1)];
        cell_list = [cell_list; local_cells];        
        
        y = 0.5*tan(new_theta);
        theta = new_theta;
        k = k+1;
    end
    % node_list = [all_x',all_y'];
    save("cell_list.mat", "cell_list");
    save("results/region_btm_cells.txt", "cell_list", "-ascii");
    save("results/region_btm_xlist.txt", "xlist", "-ascii");
    save("results/region_btm_tlist.txt", "tlist", "-ascii");
end
