function [LA_eig, LA_eigf, LA_A, LA_M, A_inner_xx, A_inner_xy, A_inner_yy, bd_dof_idx] = calc_eigen_bounds_any_order_for_quotients(tri_intval,N,ord)
    
    neig = 3;
    mesh=get_mesh_for_cg(tri_intval,N);
    vert = mesh.nodes;
    edge = mesh.edges;
    tri  = mesh.elements;
    bd   = mesh.edges(is_edge_on_bdry(mesh,N),:);
    
    Lagrange_order = ord;
    [LA_eig, LA_eigf, ~, LA_A, LA_M, A_inner_xx, A_inner_xy, A_inner_yy, bd_dof_idx] = Lagrange_dirichlet_eig(Lagrange_order, vert, edge, tri, bd, neig);

    
end
