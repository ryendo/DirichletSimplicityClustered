function [new_xlist,new_ylist,new_theta] = get_x_nodex_for_same_theta(xlist,ylist)

    a1 = xlist(1);
    b1 = ylist(1);
    theta = atan2(b1,a1);

    dt_list = [];
    N = length(xlist);

    for k = 1:(N-1)
        a1 = xlist(k);
        b1 = ylist(k);
        a3 = xlist(k+1);
        b3 = ylist(k+1);
        
        if a3^2+b3^2 > 1.00 %the last point in xy list is not trusted.
            break
        end

        dt = get_dx_dt_for_current_p1(a1,b1,a3,b3);
        dt_list = [dt_list,dt];

    end

    k = max(k,N-1);

    new_theta = theta + min(dt_list);
    new_xlist = xlist(1:k+1);
    new_ylist = xlist * tan(new_theta);

end