function [cell_list, cell_list_by_idx, xlist, tlist]=get_node_list()

    [xlist,ylist_base]=get_x_nodex_for_same_y(0.5);
    
    N = length(xlist)
    tlist = zeros(N,1);
    cell_list = []; % x1,x2,t1,t2;
    cell_list_by_idx = []; % (x1_idx,x2_idx,t1_idx,t2_idx);

    a = xlist(end); b=ylist_base(end);
    t = atan2(b,a);
    tlist(1) = t;
    t_idx = 2;
    for k = N-1:-1:1
        a2 = xlist(k);b2=ylist_base(k);
        a3 = xlist(k+1);b3=ylist_base(k+1);
        t1 = atan2(b3,a3);
        t2 = atan2(b2,a2);
        tlist(t_idx) = t2;  

        if k == 1
            display("left corner cell index:")
            size(cell_list,1)+1
        end

        for idx = k:N-1
            x1 = xlist(idx);
            x2 = xlist(idx+1);
            cell_list = [cell_list; x1,x2,t1,t2];
            cell_list_by_idx = [cell_list_by_idx; idx,idx+1,t_idx, t_idx+1];

        end
        t_idx = t_idx +1;
    end
    size(cell_list_by_idx)
    t_idx = t_idx -1;

    t = atan2(ylist_base(1), xlist(1));
    ylist = tan(t)*xlist;
    y = ylist(1);
    k=1;
    while y < 0.1
        [new_xlist,new_ylist,new_theta] = get_x_nodex_for_same_theta(xlist,ylist);
        xlist = new_xlist;
        ylist = new_ylist;
        local_cells = [xlist(1:end-1)',xlist(2:end)',t*ones(N-1,1),new_theta*ones(N-1,1)];
        cell_list = [cell_list; local_cells];
        local_regions_by_idx = [(1:N-1)',(2:N)', t_idx*ones(N-1,1), (t_idx+1)*ones(N-1,1)];
        cell_list_by_idx = [cell_list_by_idx; local_regions_by_idx];
        
        y = ylist(1);
        tlist(t_idx+1) = new_theta;
        k = k+1;
    end
    % node_list = [all_x',all_y'];
    save("results/region_btm_cells.txt", "cell_list", "-ascii");
    save("results/region_btm_cells_by_idx.txt", "cell_list_by_idx", "-ascii");
    save("results/region_btm_xlist.txt", "xlist", "-ascii");
    save("results/region_btm_tlist.txt", "tlist", "-ascii");
end