function eig_bounds = cell_lower_eig_bound(cell_data)

    % Extract geometry    
    x1 = I_intval(str2num(cell_data.x_inf));
    x2 = I_intval(str2num(cell_data.x_sup));
    t1 = I_intval(str2num(cell_data.theta_inf));
    t2 = I_intval(str2num(cell_data.theta_sup));
    
    % Calculate vertices
    a1 = x1;  b1 = x1 * tan(t1);
    a2 = x1;  b2 = x1 * tan(t2);
    a3 = x2;  b3 = x2 * tan(t1);
    a4 = x2;  b4 = x2 * tan(t2);
    
    % Extract mesh/FEM parameters    
    mesh_size_rho = cell_data.mesh_size_lower_cr;
    
    isLG = I_mid(cell_data.isLG);

    if isLG
        % =========================================================
        % Lehmann-Goerisch (LG) Method Branch
        % This computes high-accuracy lower bounds.
        % =========================================================
        
        mesh_size_LG = cell_data.mesh_size_lower_LG;
        
        % Step 1: Preliminary CR calculation to get shift parameter 'rho'
        % -------------------------------------------------------------
        % Note: Hardcoded parameters for the preliminary step as per original code
        % or utilize parameters if needed. Here we follow the snippet logic.
        neig = 3;
        
        % Geometry for mesh generation (using outer corner a4, b4)
        a = a4;
        b = b4;
        
        mesh_rho = make_mesh_by_gmsh(a, b, mesh_size_rho);
        vert_rho = mesh_rho.nodes;
        edge_rho = mesh_rho.edges;
        tri_rho  = mesh_rho.elements;
        bd_rho   = mesh_rho.boundary_edges;
        is_bnd = ismember(edge_rho, bd_rho, 'rows');

        global INTERVAL_MODE
        tri_by_edge = find_tri2edge(tri_rho, edge_rho);
        bd_edge_ids = find(is_bnd > 0);
        
        [A0, A1] = create_matrix_crouzeix_raviart(tri_rho, edge_rho, vert_rho, tri_by_edge);
        
        ne = size(edge_rho, 1);
        dof_idx = 1:ne;
        dof_idx(bd_edge_ids) = [];
        
        CR_A0 = A0(dof_idx, dof_idx);
        CR_A1 = A1(dof_idx, dof_idx);
        
        hmax = find_mesh_hmax(vert_rho, edge_rho);
        
        if INTERVAL_MODE
            CR_eig = veigs(CR_A1, CR_A0, neig + 1, 'sm');
        else
            CR_eig = eigs(CR_A1, CR_A0, neig + 1, 'sm');
        end
        
        Ch_val = I_intval('0.1893') * hmax;
        
        % CR-based lower bounds
        temp_bounds = CR_eig ./ (1 + CR_eig .* (Ch_val^2));
        [~, idx] = sort(I_mid(temp_bounds));
        temp_bounds = temp_bounds(idx);

        % Select rho (shift) from the 4th eigenvalue approximation
        rho = temp_bounds(4);
        
        % Step 2: High-order Lehmann-Goerisch Setup
        % -------------------------------------------------------------
        % Define Triangle for LG
        tri_intval = [I_intval('0'), I_intval('0'); ...
                      I_intval('1'), I_intval('0'); ...
                      a,             b];
                      
        % Make mesh for LG
        mesh_LG = make_mesh_by_gmsh(a, b, mesh_size_LG);
        vert_LG = mesh_LG.nodes;
        edge_LG = mesh_LG.edges;
        tri_LG  = mesh_LG.elements;
        bd_LG   = mesh_LG.boundary_edges;
        is_bnd_lg = ismember(edge_LG, bd_LG, 'rows');
        
        
        Lagrange_order = str2num(cell_data.fem_order_lower_LG);

        % Compute High-Order Upper Bounds & Eigenfunctions (Lagrange)
        [LA_eig, LA_eigf, LA_eigf_with_bdry, LA_A, LA_M, ~, ~, ~, ~] = ...
            laplace_eig_lagrange_detailed(Lagrange_order, vert_LG, edge_LG, tri_LG, bd_LG, neig);

        % Step 3: Construct LG Matrices
        % -------------------------------------------------------------
        A0 = LA_eigf' * LA_A * LA_eigf;
        A1 = LA_eigf' * LA_M * LA_eigf;
        
        RT_order = Lagrange_order;
        
        % Solve H(div) problem using Raviart-Thomas (RT) elements
        [mat_pih, RTRT] = RT_Hdiv_problem_dirichlet(RT_order, Lagrange_order, vert_LG, edge_LG, tri_LG, bd_LG, LA_eigf_with_bdry);

        A2 = mat_pih' * RTRT * mat_pih;
        
        % Matrices for Generalized Eigenvalue Problem
        AL = A0 - rho * A1;
        BL = A0 - 2 * rho * A1 + rho * rho * A2;
        
        % Step 4: Solve Interval Generalized Eigenvalue Problem
        % -------------------------------------------------------------
        % Symmetrize via Interval Hull
        AL = I_hull(AL, AL'); 
        BL = I_hull(BL, BL');

        mus = I_eig(AL, BL, neig);

        % Apply Lehmann-Goerisch Theorem
        LG_eig_low = rho - rho ./ (1 - mus(end:-1:1));
        
        % Sort and extract
        [~, idx] = sort(I_mid(LG_eig_low));
        LG_eig_low = LG_eig_low(idx);
        
        eig_bounds_ = I_intval(LG_eig_low);
        [~, idx] = sort(I_mid(eig_bounds_));
        eig_bounds = I_inf(eig_bounds_(idx)');

    else    
        % =========================================================
        % Standard CR Mode (isLG == 0)
        % =========================================================
        neig = 3;
        a_ = a4;
        b_ = b4;
        
        mesh_rho = make_mesh_by_gmsh(a_, b_, mesh_size_rho);
        vert_rho = mesh_rho.nodes;
        edge_rho = mesh_rho.edges;
        tri_rho  = mesh_rho.elements;
        bd_rho   = mesh_rho.boundary_edges;
        is_bnd = ismember(edge_rho, bd_rho, 'rows');

        global INTERVAL_MODE
        tri_by_edge = find_tri2edge(tri_rho, edge_rho);
        bd_edge_ids = find(is_bnd > 0);
        
        [A0, A1] = create_matrix_crouzeix_raviart(tri_rho, edge_rho, vert_rho, tri_by_edge);
        
        ne = size(edge_rho, 1);
        dof_idx = 1:ne;
        dof_idx(bd_edge_ids) = [];
        
        CR_A0 = A0(dof_idx, dof_idx);
        CR_A1 = A1(dof_idx, dof_idx);
        
        hmax = find_mesh_hmax(vert_rho, edge_rho);
        
        if INTERVAL_MODE
            CR_eig = veigs(CR_A1, CR_A0, neig + 1, 'sm');
        else
            CR_eig = eigs(CR_A1, CR_A0, neig + 1, 'sm');
        end
        
        Ch_val = I_intval('0.1893') * hmax;
        
        % Compute validated lower bounds using CR error estimate
        eig_bounds_full = CR_eig ./ (1 + CR_eig .* (Ch_val^2));

        eig_bounds = I_inf(eig_bounds_full(1:3));
    end
end