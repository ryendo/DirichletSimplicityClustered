function eig_lg = calc_eigen_bounds_any_order_1k(tri_intval,N,ord,k)
    
    
    neig = k;

    % Compute upper and lower bounds
    cr_N = 32;
    Ch=intval('0.1893')/cr_N;
    cr_llams = calc_eigen_bounds_cr(tri_intval,cr_N,1:(neig+1));
    rho = max(cr_llams); rho = rho/(1+rho*Ch^2);

    % Create mesh for L-G method
    mesh=get_mesh(tri_intval,N);
    vert = mesh.nodes;
    edge = mesh.edges;
    tri  = mesh.elements;
    bd   = mesh.edges(is_edge_on_bdry(mesh,N),:);
    
    Lagrange_order = ord;
    [LA_eig, LA_eigf, LA_eigf_with_bdry, LA_A, LA_M, A_inner_xx, A_inner_xy, A_inner_yy, bd_dof_idx, f_for_Hdiv_problem, point_for_Hdiv_problem] = Lagrange_dirichlet_eig(Lagrange_order, vert, edge, tri, bd, neig);

    disp('compute dirichlet eig lower bound by L-G thm')
    
    A0 = LA_eigf' * LA_A * LA_eigf;
    A1 = LA_eigf' * LA_M * LA_eigf;
    
    RT_order = ord;
    mat_pih = [];
    for i = 1:neig
        [pih, RTRT] = RT_Hdiv_problem_dirichlet(RT_order, Lagrange_order, vert, edge, tri, bd, LA_eigf_with_bdry(:,i));
        mat_pih = [mat_pih, pih];
    end
    A2 = mat_pih' * RTRT * mat_pih;
    
    %A2_diff = A_lg - A2
    
    AL = A0 - rho * A1;
    BL = A0 - 2*rho*A1 + rho*rho*A2;
    
    AL=(AL+AL')/2; BL=(BL+BL')/2;
    eigsALBL = veig(AL, BL,1:neig);
    [~,idx] = sort(mid(eigsALBL));
    mus = eigsALBL(idx);
    LG_eig_low = rho - rho./(1-mus(end:-1:1));
    
    % eig_cr = hull(CR_eig, CR_eig_low);
    hull_up_low = hull(LA_eig, LG_eig_low);
    eig_lg = hull_up_low(1:k);
end
