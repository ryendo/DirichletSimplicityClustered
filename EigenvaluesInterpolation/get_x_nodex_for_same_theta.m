function [new_xlist,new_theta] = get_x_nodex_for_same_theta(xlist,theta)

    dt_list = [];
    N = length(xlist);

    for k = 1:(N-1)
        a1 = xlist(k);
        b1 = a1*tan(theta);
        a3 = xlist(k+1);
        b3 = a3*tan(theta);
        
        if a3^2+b3^2 > 1.00 %the last point in xy list is not trusted.
            break
        end

        try
            dt = get_dx_dt_for_current_p1(a1,b1,a3,b3);
            dt_list = [dt_list,dt];

        catch ME
            % Handle the error
            fprintf("An error occurred: %s\n", ME.message);
        end

    end

    k = max(k,N-1);

    new_theta = theta + min(dt_list);
    dt_list
    min(dt_list)
    new_xlist = xlist(1:k+1);
end
