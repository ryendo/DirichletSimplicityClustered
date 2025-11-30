function region_list=get_node_list()

    [xlist,ylist_base]=get_x_nodex_for_same_y(0.5);
    
    N = length(xlist);
    region_list = []; % x1,x2,t1,t2;
    for k = N-1:-1:1
        a2 = xlist(k);b2=ylist_base(k);
        a3 = xlist(k+1);b3=ylist_base(k+1);
        t1 = atan2(b3,a3);
        t2 = atan2(b2,a2);

        for idx = k:N-1
            x1 = xlist(idx);
            x2 = xlist(idx+1);
            region_list = [region_list; x1,x2,t1,t2];
        end
    end

    t = atan2(ylist_base(1), xlist(1));
    ylist = tan(t)*xlist;
    y = ylist(1);
    all_x = xlist;
    all_y = ylist;
    k=1;
    while y<0.1
        [new_xlist,new_ylist,new_theta] = get_x_nodex_for_same_theta(xlist,ylist);
        xlist = new_xlist;
        ylist = new_ylist;
        % all_x = [all_x, xlist];
        % all_y = [all_y, ylist];
        local_regions = [xlist(1:end-1)',xlist(2:end)',t*ones(N-1,1),new_theta*ones(N-1,1)];
        region_list = [region_list; local_regions];

        y = ylist(1);
        t = new_theta;
        k = k+1;
    end
    % node_list = [all_x',all_y'];
end