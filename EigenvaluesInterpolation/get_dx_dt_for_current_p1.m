function delta_theta = get_dx_dt_for_current_p1(a1,b1,a3,b3)
%         (a4,b4)
%       /  |dy
% (a2,b2)  |
%   |     /(a3,b3)
% (a1,b1)/
%
%  GAP = v3(p1)-v2(p1)
%  v3(p4) - v2(p1) = sigma GAP (sigma ~ 0.25)
%  Since  v3(p4) = v3(p3) + Gy dy
%  v3(p3) + Gy dy  -v3(p1) +v3(p1) - v2(p1) = sigma GAP
%  i.e., (-Gy) dy = v3(p3) - v3(p1) + (1-sigma)GAP;

    h = 0.005;
    [v1,v2,v3,v4]=get_approximate_eigenvalue(a1,b1);
    [v1_hx,v2_hx,v3_hx,v4_hx]=get_approximate_eigenvalue(a1+h,b1);
    [v1_hy,v2_hy,v3_hy,v4_hy]=get_approximate_eigenvalue(a1,b1+h);
    G2_x = (v2_hx - v2)/h;
    G2_y = (v2_hy - v2)/h;
    G3_x = (v3_hx - v3)/h;
    G3_y = (v3_hy - v3)/h;
    v3_p1 = v3;

    GAP = v3-v2;
    
    min_with_ratio = 0.25;
    theta = atan2(b1,a1);

    [v1_p3,v2_p3,v3_p3,v4_p3]=get_approximate_eigenvalue(a3,b3);
    delta_y_initial = - ((v3_p3 - v3_p1)+(1-min_with_ratio)*GAP) / G3_y;

    f = @(t) predict_diag_gap_for_p1(t, a1, b1, a3, b3, min_with_ratio);
    delta_y = fzero(f, [delta_y_initial/5 delta_y_initial*10]);
    
    a4 = a3;
    b4 = b3+delta_y;

    theta2 = atan2(b4,a4);
    delta_theta = theta2 - theta;
    [v1_p1, v2_p1, v3_p1, v4_p1] = get_approximate_eigenvalue(a1,b1);
    [v1_p4, v2_p4, v3_p4, v4_p4] = get_approximate_eigenvalue(a4,b4);

    GAP = v3_p1 -v2_p1;
    v23_gap_validation = v3_p4 - v2_p1;

    if v23_gap_validation > min_with_ratio*GAP*0.85 
        rel_width=v23_gap_validation/GAP;
        rel_width;
        % display("OK");
    else
        [v23_gap_validation,GAP/2];
        rel_width=v23_gap_validation/GAP;
        rel_width
        display('NG')
    end
    %[a1,b1,a4,b4]

end
