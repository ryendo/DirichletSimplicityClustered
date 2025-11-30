function gap_diff_from_target = predict_diag_gap_for_p1(delta_y,a1,b1,a3,b3,target_ratio)
    [v1,v2,v3,v4] = get_approximate_eigenvalue(a1,b1);

    GAP0 = v3-v2;

    a4 = a3;
    b4 = b3+delta_y;
    v2_p1 = v2;
    [v1_p4, v2_p4, v3_p4, v4_p4] = get_approximate_eigenvalue(a4,b4);

    new_gap = v3_p4 - v2_p1;
    gap_diff_from_target = new_gap/GAP0 - target_ratio;
    
end
