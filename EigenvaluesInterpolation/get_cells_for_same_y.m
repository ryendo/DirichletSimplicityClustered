function cell_list = get_cells_for_same_y(x_start)
    y0=tan(pi/60)/2-(1E-4);
    x = x_start;
    y = y0;
    theta2 = atan2(y,x);

    xlist = []; ylist=[];
    cell_list = [];

    % (x,y0)  --------- (x + dx)
    %  |              / |   / 
    %  |            /   (next_p1_x,next_p1_y)
    %  (p1_x,p1_y)

    while x <= 1

        [dx, ~] = get_dx_dt_for_current_p2(x,y0);
        dx = 0.75 * dx;
        next_p1_x = x+dx;
        if next_p1_x > 1
            next_p1_x = 1.001;
        end

        theta1 = atan2(y0, next_p1_x);
        p1_y = x*tan(theta1);

        cell_list = [cell_list; x, next_p1_x, theta1, theta2];

        x = next_p1_x;
        theta2 = theta1; % next thetat2 is current theta1

    end

end
