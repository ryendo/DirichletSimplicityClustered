function [xlist,ylist] = get_x_nodes_for_same_y(x_start)
    a0=x_start;
    b0=tan(pi/60)/2-(1E-4);
    x = x_start;
    [dx,dt] = get_dx_dt_for_current_p2(x,b0);

    x2 = x+dx;
    theta1 = atan2(b0,x2);
    b1 = x*tan(theta1);
    xlist = [x];
    ylist = [b1];

    while x <= 1
        x = x+dx;
	if x > 1.0001
           x = 1.0001;
           xlist = [xlist,x];
           ylist = [ylist,b0];
	   break
	end
        xlist = [xlist,x];
        [dx,dt] = get_dx_dt_for_current_p2(x,b0);
        theta1 = atan2(b0,x+dx);
        b1 = x*tan(theta1);
        ylist = [ylist,b1];
    end
end
