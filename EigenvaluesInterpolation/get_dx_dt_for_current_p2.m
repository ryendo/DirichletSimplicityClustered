function [delta_x,delta_theta] = get_dx_dt_for_current_p2(a2,b2)
    h = 0.005;
    [v1,v2,v3,v4]=get_approximate_eigenvalue(a2,b2);
    [v1_hx,v2_hx,v3_hx,v4_hx]=get_approximate_eigenvalue(a2+h,b2);
    [v1_hy,v2_hy,v3_hy,v4_hy]=get_approximate_eigenvalue(a2,b2+h);
    G2_x = (v2_hx - v2)/h;
    G2_y = (v2_hy - v2)/h;
    G3_x = (v3_hx - v3)/h;
    G3_y = (v3_hy - v3)/h;

    GAP = v3-v2;
    
    min_with_ratio = 0.25;
    theta2 = atan2(b2,a2);

    delta_x_initial = - (1-min_with_ratio)*GAP / (G2_y*tan(theta2)+G3_x + G3_y*tan(theta2));

    f = @(t) predict_diag_gap_for_p2(t, a2, b2, min_with_ratio);
    delta_x = fzero(f, [delta_x_initial/5 delta_x_initial*10]);

    theta1 = atan2(b2,a2+delta_x);
    delta_theta = theta2 - theta1;
    a1 = a2;
    b1 = a1*tan(theta1);
    a3 = a2+delta_x;
    b3 = a3*tan(theta1);
    a4 = a3;
    b4 = a4*tan(theta2);
    [v1_p1, v2_p1, v3_p1, v4_p1] = get_approximate_eigenvalue(a1,b1);
    [v1_p4, v2_p4, v3_p4, v4_p4] = get_approximate_eigenvalue(a4,b4);

    GAP = v3_p1 -v2_p1;
    v23_gap_validation = v3_p4 - v2_p1;

    if v23_gap_validation > min_with_ratio*GAP*0.85 
        rel_width=v23_gap_validation/GAP;
        rel_width
        display("OK")
    else
        [v23_gap_validation,GAP/2];
        rel_width=v23_gap_validation/GAP;
        rel_width
        display('NG')
    end
    [a1,b1,a4,b4]

end