function eig_bounds = calc_eigen_bounds_any_order_lam1(tri_intval,N_LG,N_rho,ord,isLG)

    if isLG==1
    
        neig = 1;
    
        % Compute upper and lower bounds
        
        mesh_rho=get_mesh_for_cg(tri_intval,N_rho);
        vert_rho = mesh_rho.nodes;
        edge_rho = mesh_rho.edges;
        tri_rho  = mesh_rho.elements;
        bd_rho   = mesh_rho.edges(is_edge_on_bdry(mesh_rho,N_rho),:);
        Ch=I_intval('0.493')/N_rho;
        cg_llams = Lagrange_dirichlet_eig_for_rho(1, vert_rho, edge_rho, tri_rho, bd_rho, neig+1);
        % cr_llams = calc_eigen_bounds_cr(tri_intval,cr_N,1:(neig+1));
        rho = max(cg_llams); rho = rho/(1+rho*Ch^2);
        
        
        % Create mesh for L-G method
        mesh=get_mesh(tri_intval,N_LG);
        vert = mesh.nodes;
        edge = mesh.edges;
        tri  = mesh.elements;
        bd   = mesh.edges(is_edge_on_bdry(mesh,N_LG),:);
        
        Lagrange_order = ord;
        [LA_eig, LA_eigf, LA_eigf_with_bdry, LA_A, LA_M, A_inner_xx, A_inner_xy, A_inner_yy, bd_dof_idx] = Lagrange_dirichlet_eig(Lagrange_order, vert, edge, tri, bd, neig);
        
    
        % disp('compute dirichlet eig lower bound by L-G thm')
        
        A0 = LA_eigf' * LA_A * LA_eigf;
        A1 = LA_eigf' * LA_M * LA_eigf;
        
        RT_order = ord;
        
        [mat_pih, RTRT] = RT_Hdiv_problem_dirichlet(RT_order, Lagrange_order, vert, edge, tri, bd, LA_eigf_with_bdry);
        A2 = mat_pih' * RTRT * mat_pih;
        
        
        %A2_diff = A_lg - A2
        
        
        AL = A0 - rho * A1;
        BL = A0 - 2*rho*A1 + rho*rho*A2;
        
        AL=(AL+AL')/2; BL=(BL+BL')/2;
        eigsALBL = I_veig(AL, BL,1:neig);
        [~,idx] = sort(I_mid(eigsALBL));
        mus = eigsALBL(idx);
        LG_eig_low = rho - rho./(1-mus(end:-1:1));
        
        % eig_cr = hull(CR_eig, CR_eig_low);
        hull_up_low = hull(intval(LA_eig), intval(LG_eig_low));
        eig_bounds = hull_up_low(1);

    else
        neig = 1;
    
        % Compute upper and lower bounds
        mesh_rho=get_mesh_for_cg_y_reduced(tri_intval,N_rho);
        vert_rho = mesh_rho.nodes;
        edge_rho = mesh_rho.edges;
        tri_rho  = mesh_rho.elements;
        bd_rho   = mesh_rho.edges(is_edge_on_bdry(mesh_rho,N_rho),:);
        Ch=I_intval('0.493')/N_rho;
        cg_lams = Lagrange_dirichlet_eig_for_rho(1, vert_rho, edge_rho, tri_rho, bd_rho, neig+1);
        % cr_llams = calc_eigen_bounds_cr(tri_intval,cr_N,1:(neig+1));
        llams = cg_lams./(1+cg_lams*Ch^2);
        hull_up_low = hull(cg_lams,llams);
        eig_bounds = hull_up_low(1);
    end
    
end
