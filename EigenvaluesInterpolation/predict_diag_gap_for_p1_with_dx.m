function gap_diff_from_target = predict_diag_gap_for_p1_with_dx(delta_x,a1,b1,target_ratio)
    [v1,v2_p1,v3_p1,v4] = get_approximate_eigenvalue(a1,b1);

    GAP0 = v3_p1-v2_p1;

    theta = atan2(b1,a1);
    a4 = a1 + delta_x;

    delta_y = delta_x; % assume dx=dy
    b4 = b1 + delta_x*tan(theta) + delta_y;

    [v1_p4, v2_p4, v3_p4, v4_p4] = get_approximate_eigenvalue(a4,b4);
    v3_p4 = v3_p4*(1-0.005);

    new_gap = v3_p4 - v2_p1;
    gap_diff_from_target = new_gap/GAP0 - target_ratio;
    
end
