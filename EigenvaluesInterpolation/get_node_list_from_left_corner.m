function cell_list=get_node_list_from_left_corner(a1,b1)

    %from left to right

    x = a1;
    y = b1;

    k=1;
    t = atan2(y,x);
    cell_list =[];
    while y < 0.1
        [k,y]
        max_x = cos(t);

        x=0.5;
        xlist = [x];
        dt_list = [];
        while x < max_x
            p1_x = x;
            p1_y = p1_x * tan(t)
            [dx,dt] = get_dx_dt_for_current_p1(p1_x, p1_y);
            x = x + dx*0.9;
            xlist= [xlist; x];
            dt_list = [dt_list; dt];
        end
        dt = min(dt_list);

        nt = t + dt;

        N = size(xlist,1)-1;
        local_cells = [xlist(1:end-1),xlist(2:end),t*ones(N,1),nt*ones(N,1)];
        cell_list = [cell_list; local_cells];        
        
        t = nt;
        y = 0.5*tan(t);
        k = k+1;
    end
    save("cell_list_from_y_0.026.mat", "cell_list");
end
