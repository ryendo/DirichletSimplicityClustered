function cell_ub = cell_upper_eig_bound(cell_data)
    % Extract and convert geometry parameters to intervals
    x1 = I_intval(str2num(cell_data.x_inf));
    x2 = I_intval(str2num(cell_data.x_sup));
    t1 = I_intval(str2num(cell_data.theta_inf));
    t2 = I_intval(str2num(cell_data.theta_sup));
    
    % Define vertices
    % Note: a1, b1 corresponds to x1, t1
    a1 = x1;  b1 = x1 * tan(t1);
    
    % Setup FEM parameters
    neig = 3;
    fem_ord = I_mid(I_intval(cell_data.fem_order_upper));
    mesh_size = cell_data.mesh_size_upper;
    
    % Mesh Generation (using Gmsh)
    mesh_rho = make_mesh_by_gmsh(a1, b1, mesh_size);        
    vert_rho = I_intval(mesh_rho.nodes);
    edge_rho = mesh_rho.edges;
    tri_rho  = mesh_rho.elements;
    bd_rho   = mesh_rho.boundary_edges;
    
    % Compute upper bounds using Lagrange FEM
    cg_lams = upper_eig_bound(fem_ord, vert_rho, edge_rho, tri_rho, bd_rho, neig);

    % Return the supremum of the first 3 eigenvalues
    cell_ub = I_sup(cg_lams(1:3));
end