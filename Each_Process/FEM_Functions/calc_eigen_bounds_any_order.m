% calc_eigen_bounds_any_order: Computes rigorous upper and lower bounds for the eigenvalues of the Dirichlet Laplacian.
%
% This function can operate in two modes, controlled by the 'isLG' flag.
% 1. Lehmann-Goerisch (LG) mode (isLG=1): Provides high-accuracy lower bounds using the Lehmann-Goerisch theorem.
% 2. Standard CG mode (isLG~=1): Computes upper bounds via the continuous Galerkin (CG) finite element method 
%    and lower bounds using the error estimate.
%
% All calculations are performed using interval arithmetic to ensure mathematical rigor.


function eig_bounds = calc_eigen_bounds_any_order(tri_intval,N_LG,N_rho,ord,isLG)

    if isLG==1
        % --- Lehmann-Goerisch (LG) Method Branch ---
        % This branch computes high-accuracy lower bounds.
        
        % Define the number of eigenvalues to compute.
        neig = 3;
    
        % --- Step 1: Compute the constant 'rho'. ---
        % This is required by the Lehmann-Goerisch theorem.

        % Generate a right-triangular mesh for the Continuous Galerkin (CG) method.
        mesh_rho=get_mesh_for_cg(tri_intval,N_rho);
        vert_rho = mesh_rho.nodes;
        edge_rho = mesh_rho.edges;
        tri_rho  = mesh_rho.elements;
        bd_rho   = mesh_rho.edges(is_edge_on_bdry(mesh_rho,N_rho),:);

        % Define the constant for the error estimate.
        Ch=I_intval('0.493')/N_rho;

        % Compute rough upper bounds for the eigenvalues using the P1 Lagrange FEM.
        cg_llams = Lagrange_dirichlet_eig_for_rho_vectorized(1, vert_rho, edge_rho, tri_rho, bd_rho, neig+1);
        rho = max(cg_llams); rho = rho/(1+rho*Ch^2);
        
        % --- Step 2: Set up the main high-order Lehmann-Goerisch method. ---
        % Create a mesh for the LG computation.
        mesh=get_mesh(tri_intval,N_LG);
        vert = mesh.nodes;
        edge = mesh.edges;
        tri  = mesh.elements;
        bd   = mesh.edges(is_edge_on_bdry(mesh,N_LG),:);
        
        % Set the polynomial order for the Lagrange elements.
        Lagrange_order = ord;

        % Compute high-accuracy upper bounds (LA_eig) and corresponding eigenfunctions (LA_eigf)
        % using the high-order Lagrange method. Also get the system matrices A and M.
        
        [LA_eig, LA_eigf, LA_eigf_with_bdry, LA_A, LA_M, ~, ~, ~, ~] = Lagrange_dirichlet_eig_vectorized(Lagrange_order, vert, edge, tri, bd, neig);

        % --- Step 3: Construct the matrices for the Lehmann-Goerisch generalized eigenvalue problem. ---
        % The problem is of the form AL*x = mu*BL*x.
        
        % Form the matrices A0 and A1 from the high-order system matrices and eigenfunctions.
        A0 = LA_eigf' * LA_A * LA_eigf;
        A1 = LA_eigf' * LA_M * LA_eigf;
        
        % Set the order for the Raviart-Thomas (RT) elements.
        RT_order = ord;
        
        % Solve the H(div) problem using RT elements to construct the matrix A2.
        % This involves projecting the Lagrange eigenfunctions onto the RT space.
        [mat_pih, RTRT] = RT_Hdiv_problem_dirichlet(RT_order, Lagrange_order, vert, edge, tri, bd, LA_eigf_with_bdry);
        A2 = mat_pih' * RTRT * mat_pih;
        
        % Construct the matrices AL and BL for the generalized eigenvalue problem.
        AL = A0 - rho * A1;
        BL = A0 - 2*rho*A1 + rho*rho*A2;
        
        % --- Step 4: Solve the generalized eigenvalue problem and compute the final bounds. ---
        % Ensure matrices are symmetric using the interval hull operation.
        AL=I_hull(AL,AL'); BL=I_hull(BL,BL');

        % Solve the interval generalized eigenvalue problem AL*x = mu*BL*x for eigenvalues 'mu'.
        mus =  I_veig(AL, BL,1:neig);

        % Calculate the rigorous lower bounds for the original problem's eigenvalues
        % using the 'mus' and the reference 'rho'. The formula is derived from the LG theorem.
        LG_eig_low = rho - rho./(1-mus(end:-1:1));
        
        % Combine the high-accuracy upper bounds (LA_eig) and the new lower bounds (LG_eig_low)
        % into a final interval vector.
        [~,idx] = sort(I_mid(LG_eig_low));
        LG_eig_low = LG_eig_low(idx);
        [~,idx] = sort(I_mid(LA_eig));
        LA_eig = LG_eig_low(idx);

        hull_up_low = I_hull(I_intval(LA_eig), I_intval(LG_eig_low));

        % Extract the desired eigenvalue bounds.
        eig_bounds_ = hull_up_low(2:3);

        % Sort the computed eigenvalues.
        [~,idx] = sort(I_mid(eig_bounds_));
        eig_bounds = eig_bounds_(idx);

    else
        % --- Standard CG Method Branch ---
        % This branch computes bounds using a standard a posteriori error estimate,
        % which is faster but provides less accurate lower bounds than the LG method.

        % Define the number of eigenvalues to compute.
        neig = 3;
    
        % --- Step 1: Compute upper and lower bounds. ---
        % Generate a mesh suitable for the CG method.
        if tri_intval(6)<0.09
            N_v = 10;
            mesh_rho=get_mesh_for_cg_y_reduced(tri_intval,N_rho,N_v);
        else
            mesh_rho=get_mesh_for_cg(tri_intval,N_rho);
        end
        
        vert_rho = mesh_rho.nodes;
        edge_rho = mesh_rho.edges;
        tri_rho  = mesh_rho.elements;
        bd_rho   = mesh_rho.edges(is_edge_on_bdry(mesh_rho,N_rho),:);

        % Define the constant for the a posteriori error estimate.
        Ch=I_intval('0.493')/N_rho;

        % Compute the upper bounds for the eigenvalues using P1 Lagrange FEM.
        cg_lams = Lagrange_dirichlet_eig_for_rho_vectorized(1, vert_rho, edge_rho, tri_rho, bd_rho, neig+1);
        
        % Compute the corresponding lower bounds using the error formula:
        llams = cg_lams./(1+cg_lams*Ch^2);

        % Combine the upper bounds (cg_lams) and lower bounds (llams) into an interval vector.
        hull_up_low = I_hull(cg_lams,llams);

        % Extract the desired eigenvalue bounds.
        eig_bounds = hull_up_low(2:3);
    end
    
end
