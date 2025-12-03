function [delta_x, delta_theta] = get_dx_dt_for_current_p1(a1,b1)
%         (a4,b4)
%       /  |dy
% (a2,b2)  |
%   |     /(a3,b3)
% (a1,b1)/  a3=a1+dx
%
%  GAP = v3(p1)-v2(p1)
%  LB(v3(p4)) - v2(p1) = sigma GAP  (here, sigma ~ 0.25)
%  Since  v3(p4) = v3(p1) +Gx dx + Gy dy =  v3(p1) + (Gx+Gy)*dx (here, assume dx=dy)
%  LB(v3(p1) + (Gx+Gy)dx) - v2(p1) = sigma GAP
%  Here, LB(lambda) = lambda*(1-0.005), here assume the lower bound gap is less than 0.5%.
%  i.e., (1-0.005)(v3(p1) + (Gx+Gy)*dx) - v2(p1) = sigma GAP;
%  i.e., v3(p1) + (Gx+Gy)*dx - 0.005(v3(p1) + (Gx+Gy)*dx) - v2(p1) = sigma GAP;
%  i.e., GAP + 0.995(Gx+Gy)*dx - 0.005*v3(p1)  = sigma GAP;
%  i.e.,  (Gx+Gy)dx = ( (sigma-1)GAP+0.005*v3(p1) )/0.995 ;
%  i.e.,  dx =   ( (sigma-1)GAP+0.005*v3(p1) )/0.995 / (Gx+Gy);

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
    
    min_with_ratio = 0.2; %sigma
    theta = atan2(b1,a1);

    delta_x_initial = ((min_with_ratio-1)*GAP + 0.005*v3)/0.995/(G3_y+G3_x);

    try
        f = @(t) predict_diag_gap_for_p1_with_dx(t, a1, b1, min_with_ratio);
        delta_x = fzero(f, [delta_x_initial/10 delta_x_initial*10]);
    catch ME
        disp("failed to find optimal delta_y.")
        fprintf("An error occurred: %s\n", ME.message);
    end
    delta_theta = atan2(b1 + delta_x*tan(theta) + delta_x, a1 + delta_x) - theta;

end
