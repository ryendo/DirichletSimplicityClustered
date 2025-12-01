% calc_eigen_bounds_any_order: Computes rigorous upper and lower bounds for the eigenvalues of the Dirichlet Laplacian.
%
% This function can operate in two modes, controlled by the 'isLG' flag.
% 1. Lehmann-Goerisch (LG) mode (isLG=1): Provides high-accuracy lower bounds using the Lehmann-Goerisch theorem.
% 2. Standard CG mode (isLG~=1): Computes upper bounds via the continuous Galerkin (CG) finite element method 
%    and lower bounds using the error estimate.
%
% All calculations are performed using interval arithmetic to ensure mathematical rigor.


function eig_bounds = cell_lower_eig_bound(region_cell)

        x1 = region_cell(1);
        x2 = region_cell(2);
        t1 = region_cell(3);
        t2 = region_cell(4);
        a1 = x1;  b1 = x1*tan(t1);
        a2 = x1;  b2 = x1*tan(t2);
        a3 = x2;  b3 = x2*tan(t1);
        a4 = x2;  b4 = x2*tan(t2);

        % --- Standard CG Method Branch ---
        % This branch computes bounds using a standard a posteriori error estimate,
        % which is faster but provides less accurate lower bounds than the LG method.

        % Define the number of eigenvalues to compute.
        neig = 3;
    
        % --- Step 1: Compute upper and lower bounds. ---
        % Generate a mesh suitable for the CG method.
        N_v = 5;
        % mesh_rho=get_mesh_for_cg_y_reduced(tri_intval,N_rho,N_v);
        N_v = 5;
        a_ = a4;
        b_ = b4;
        mesh_rho = make_mesh_by_gmsh(a_, b_, 0.002);
        vert_rho = mesh_rho.nodes;
        edge_rho = mesh_rho.edges;
        tri_rho  = mesh_rho.elements;
        bd_rho   = mesh_rho.boundary_edges;
        is_bnd = ismember(edge_rho, bd_rho, 'rows');


        is_bnd = ismember(edge_rho, bd_rho, 'rows');

	vert_rho = I_intval(vert_rho);

        % Define the constant for the a posteriori error estimate.

        disp('--- Compute Laplacian eigenvalues using CR element ---');

        global INTERVAL_MODE
        tri_by_edge = find_tri2edge(tri_rho, edge_rho);
        bd_edge_ids = find(is_bnd>0);
        [A0, A1] = create_matrix_crouzeix_raviart(tri_rho, edge_rho, vert_rho, tri_by_edge);
        ne = size(edge_rho,1);
        dof_idx = 1:ne;
        dof_idx(bd_edge_ids) = [];
        CR_A0 = A0(dof_idx,dof_idx);
        CR_A1 = A1(dof_idx,dof_idx);
        hmax = find_mesh_hmax(vert_rho,edge_rho);
        if INTERVAL_MODE
        CR_eig = veigs(CR_A1, CR_A0, neig+1, 'sm');
        else
        CR_eig = eigs(CR_A1, CR_A0, neig+1, 'sm');
        end
        Ch_val = I_intval(0.1893)*hmax;
        
        disp('Compute validated lower bounds (CR-based theorem)');
        eig_bounds = CR_eig ./ (1 + CR_eig .* (Ch_val^2));

        eig_bounds = eig_bounds(1:4);


    end
   
