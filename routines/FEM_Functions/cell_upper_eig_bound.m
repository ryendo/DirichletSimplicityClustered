% calc_eigen_bounds_any_order: Computes rigorous upper and lower bounds for the eigenvalues of the Dirichlet Laplacian.
%
% This function can operate in two modes, controlled by the 'isLG' flag.
% 1. Lehmann-Goerisch (LG) mode (isLG=1): Provides high-accuracy lower bounds using the Lehmann-Goerisch theorem.
% 2. Standard CG mode (isLG~=1): Computes upper bounds via the continuous Galerkin (CG) finite element method 
%    and lower bounds using the error estimate.
%
% All calculations are performed using interval arithmetic to ensure mathematical rigor.


function eig_bounds = cell_upper_eig_bound(region_cell,mesh_h)

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

        a_ = a1;
        b_ = b1;
        mesh_rho = make_mesh_by_gmsh(a_, b_, mesh_h);
        vert_rho = mesh_rho.nodes;
        edge_rho = mesh_rho.edges;
        tri_rho  = mesh_rho.elements;
        bd_rho   = mesh_rho.boundary_edges;
        is_bnd = ismember(edge_rho, bd_rho, 'rows');

        disp("Element number:")
        disp(size(tri_rho,1))
        disp("Node number:")
        disp(size(vert_rho,1))

        vert_rho = I_intval(vert_rho);
            
        % Compute the upper bounds for the eigenvalues using P1 Lagrange FEM.
        cg_eig_upper_bound = Lagrange_upper_eig_bound(1, vert_rho, edge_rho, tri_rho, bd_rho, neig+1);

        % Extract the desired eigenvalue bounds.
        eig_bounds = I_sup(cg_eig_upper_bound(1:4));
    end
   
