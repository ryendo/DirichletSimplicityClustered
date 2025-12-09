function [LA_eig, LA_eigf, LA_A, LA_M, A_inner_xx, A_inner_xy, A_inner_yy, bd_dof_idx] = calc_eigen_bounds_any_order_for_quotients(tri_intval,ord)

    neig = 3;
    mesh_size = 0.03125;
    mesh = make_mesh_by_gmsh(a1, b1, mesh_size);        
    vert = I_intval(mesh_rho.nodes);
    edge = mesh_rho.edges;
    tri  = mesh_rho.elements;
    bd   = mesh_rho.boundary_edges;
    
    Lagrange_order = ord;
    [LA_eig, LA_eigf, ~, LA_A, LA_M, A_inner_xx, A_inner_xy, A_inner_yy, bd_dof_idx] = laplace_eig_lagrange_detailed(Lagrange_order, vert, edge, tri, bd, neig);

    
end
