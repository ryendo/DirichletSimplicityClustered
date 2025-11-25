function [eig_value, eig_func, A, M] = laplace_eig_lagrange(lagrange_order, mesh, neig)

    vert = I_intval(mesh.nodes);
    edge = mesh.edges;
    tri  = mesh.elements;
    bd_edge_ids  = mesh.bd_edge_ids;
    bd_edges = edge(bd_edge_ids,:);

    ne = size(edge, 1);
    nt = size(tri,  1);
    nb = size(bd_edges,   1);
    nv = size(vert, 1);
    if nt > 1000
        disp('Processing the mesh...');
    end
    

    [tri, tri2edge] = find_tri2edge(edge, tri, ne, nt);

    %[Liu] Get boundary DOF idx of Lagrange FEM space.
    % Note that the DOF on an edge includes the DOFs on two ends of an edge, 
    % and the DOFs inside the edge when lagrange_order > 1.
    bd_dof_list_interior = zeros(1,length(bd_edge_ids)*(lagrange_order-1));
    bd_dof_list_two_ends = zeros(1,length(bd_edge_ids)*2);
    for idx = 1:length(bd_edge_ids)
        edge_idx = bd_edge_ids(idx);
        idx_inner = (idx-1)*(lagrange_order-1)+1:idx*(lagrange_order-1); 
        edge_interiori = (edge_idx-1)*(lagrange_order-1)+1:edge_idx*(lagrange_order-1);
        bd_dof_list_interior(idx_inner) = nv + edge_interiori;
        bd_dof_list_two_ends ( (idx-1)*2+1:2*idx ) = edge(edge_idx,:); 
    end
    bd_dof_idx = [unique(bd_dof_list_two_ends),bd_dof_list_interior];

    if nt > 1000
        disp('Create matrix for CG FEM applied to Laplace ...');
    end
    [A, M] = Lagrange_laplace_eig_matrix(lagrange_order, vert, edge, tri, tri2edge);

    dof_non_bd = 1:size(A,1);
    dof_non_bd(bd_dof_idx) = [];
    A0 = A(dof_non_bd,dof_non_bd);
    M0 = M(dof_non_bd,dof_non_bd);
    
    if nt > 1000
        disp('Start to solve eigenvalue problem. A0 = lambda M0');
    end

    neig = min(neig, size(A0, 1));
    [V, D] = eigs(I_mid(A0), I_mid(M0), neig, 'sm');
    [eig_value, idx] = sort(diag(D));
    eig_func_interior = V(:, idx);
    n_dof_num = size(A,1);
    eig_func = zeros(n_dof_num,neig);
    eig_func (dof_non_bd,:) = eig_func_interior;
    
    A = full(A);
    M = full(M);

end


function [A, M] = Lagrange_laplace_eig_matrix(lagrange_order, vert, edge, tri, tri2edge)
    [basis, nbasis] = Lagrange_basis(lagrange_order);
    M_ip_elem  = Lagrange_inner_product_L1L2L3_all(lagrange_order);
    M_ip_edge1 = Lagrange_inner_product_edge_L1L2L3_all(lagrange_order, 1);
    M_ip_edge2 = Lagrange_inner_product_edge_L1L2L3_all(lagrange_order, 2);
    M_ip_edge3 = Lagrange_inner_product_edge_L1L2L3_all(lagrange_order, 3);
    
    M_ip_basis_grad_ijT_all = cell(nbasis);
    M_ip_basis_ij = I_zeros(nbasis, nbasis);
    M_ip_basis_edge_ij_all = cell(3, 1);
    M_ip_basis_edge1_ij = I_zeros(nbasis, nbasis);
    M_ip_basis_edge2_ij = I_zeros(nbasis, nbasis);
    M_ip_basis_edge3_ij = I_zeros(nbasis, nbasis);
    for i = 1:nbasis
        for j = i:nbasis
            ei = Lagrange_create_coord_basis_grad(basis, i, lagrange_order);
            ej = Lagrange_create_coord_basis_grad(basis, j, lagrange_order);
            M_ip_basis_grad_ijT_all{j, i} = ei' * M_ip_elem * ej;
            
            ei = Lagrange_create_coord_basis(basis, i, lagrange_order);
            ej = Lagrange_create_coord_basis(basis, j, lagrange_order);
            M_ip_basis_ij(j, i) = ei' * M_ip_elem * ej;
            M_ip_basis_edge1_ij(j, i) = ei' * M_ip_edge1 * ej;
            M_ip_basis_edge2_ij(j, i) = ei' * M_ip_edge2 * ej;
            M_ip_basis_edge3_ij(j, i) = ei' * M_ip_edge3 * ej;
        end
    end
    M_ip_basis_edge_ij_all{1, 1} = M_ip_basis_edge1_ij + tril(M_ip_basis_edge1_ij, -1)';
    M_ip_basis_edge_ij_all{2, 1} = M_ip_basis_edge2_ij + tril(M_ip_basis_edge2_ij, -1)';
    M_ip_basis_edge_ij_all{3, 1} = M_ip_basis_edge3_ij + tril(M_ip_basis_edge3_ij, -1)';
    
    nv = size(vert, 1);
    ne = size(edge, 1);
    nt = size(tri,  1);
    
    % [Liu] Total dof of Lagrange FEM on a mesh is counted by the number of
    % DOFsat vertices, edges and the ones inside each element.
    ndof = nv + (lagrange_order-1)*ne + (lagrange_order-1)*(lagrange_order-2)/2*nt;
    A = I_sparse(ndof, ndof);
    M = I_sparse(ndof, ndof);
    
    for k = 1:nt
        vert_idx = tri(k,:);
        x1 = vert(vert_idx(1), 1); y1 = vert(vert_idx(1), 2);
        x2 = vert(vert_idx(2), 1); y2 = vert(vert_idx(2), 2);
        x3 = vert(vert_idx(3), 1); y3 = vert(vert_idx(3), 2);
        B = [x2-x1, x3-x1; y2-y1, y3-y1];
        Binv = [y3-y1, x1-x3; y1-y2, x2-x1] / det(B);
        A_grad = I_zeros(nbasis, nbasis);
        for i = 1:nbasis
            for j = 1:i
                A_grad(i, j) = trace(Binv' * M_ip_basis_grad_ijT_all{i, j} * Binv) * det(B);
            end
        end
        A_grad = A_grad + tril(A_grad, -1)';
        
        A_L2 = (M_ip_basis_ij + tril(M_ip_basis_ij, -1)') * det(B);
        
        edge_idx = tri2edge(k, :);
        edge_vec = vert(vert_idx([2,3,1]), :) - vert(vert_idx([3,1,2]), :);
        local_edge_start_vert = [2, 3, 1];
        map_dof_idx_l2g = zeros(nbasis, 1);
        map_dof_idx_l2g(1:3) = vert_idx;
        for i = 1:3
            if edge(edge_idx(i), 1) == vert_idx(local_edge_start_vert(i))
                map_dof_idx_l2g(3+(lagrange_order-1)*(i-1)+1:3+(lagrange_order-1)*i) = nv+(lagrange_order-1)*(edge_idx(i)-1)+1:1:nv+(lagrange_order-1)*edge_idx(i);
            else
                map_dof_idx_l2g(3+(lagrange_order-1)*(i-1)+1:3+(lagrange_order-1)*i) = nv+(lagrange_order-1)*edge_idx(i):-1:nv+(lagrange_order-1)*(edge_idx(i)-1)+1;
            end
        end
        map_dof_idx_l2g(3+(lagrange_order-1)*3+1:end) = nv + (lagrange_order-1)*ne + ...
                                                        ((lagrange_order-1)*(lagrange_order-2)/2*(k-1)+1:(lagrange_order-1)*(lagrange_order-2)/2*k);
    
        A(map_dof_idx_l2g, map_dof_idx_l2g) = A(map_dof_idx_l2g, map_dof_idx_l2g) + A_grad;
    
        M(map_dof_idx_l2g, map_dof_idx_l2g) = M(map_dof_idx_l2g, map_dof_idx_l2g) + A_L2;
    
    end
end

function M_ip = Lagrange_inner_product_L1L2L3_all(lagrange_order)
    ijk = create_ijk(lagrange_order);
    len = size(ijk, 1);
    M_ip = I_zeros(len, len);
    for p = 1:len
        for q = p:len
            pi = ijk(p, 1);
            pj = ijk(p, 2);
            pk = ijk(p, 3);
            qi = ijk(q, 1);
            qj = ijk(q, 2);
            qk = ijk(q, 3);
            M_ip(p, q) = Lagrange_integral_L1L2L3_ijk(pi+qi, pj+qj, pk+qk);
        end
    end
    M_ip = M_ip + triu(M_ip, 1)';
end

function M_ip = Lagrange_inner_product_edge_L1L2L3_all(lagrange_order, edge_idx)
    ijk = create_ijk(lagrange_order);
    len = size(ijk, 1);
    M_ip = I_zeros(len, len);
    for p = 1:len
        for q = p:len
            pi = ijk(p, 1);
            pj = ijk(p, 2);
            pk = ijk(p, 3);
            qi = ijk(q, 1);
            qj = ijk(q, 2);
            qk = ijk(q, 3);
            M_ip(p, q) = Lagrange_integral_edge_L1L2L3_ijk(pi+qi, pj+qj, pk+qk, edge_idx);
        end
    end
    M_ip = M_ip + triu(M_ip, 1)';
end

function e = Lagrange_create_coord_basis(basis, idx, Lgrange_order)
    len = Lagrange_get_nbasis(Lgrange_order);
    e = zeros(len, 1);
    
    i = basis(idx, 1);
    j = basis(idx, 2);
    k = basis(idx, 3);
    e(map_ijk_to_idx(i, j, k, Lgrange_order), 1) = 1;
    end
    
    function e = Lagrange_create_coord_basis_grad(basis, idx, Lgrange_order)
    len = Lagrange_get_nbasis(Lgrange_order);
    dudx = zeros(len, 1);
    dudy = zeros(len, 1);
    
    i = basis(idx, 1);
    j = basis(idx, 2);
    k = basis(idx, 3);
    
    dudx(map_ijk_to_idx(i,   j,   k,   Lgrange_order), 1) = -i + j;
    dudx(map_ijk_to_idx(i-1, j+1, k,   Lgrange_order), 1) = -i;
    dudx(map_ijk_to_idx(i-1, j,   k+1, Lgrange_order), 1) = -i;
    dudx(map_ijk_to_idx(i+1, j-1, k,   Lgrange_order), 1) =  j;
    dudx(map_ijk_to_idx(i,   j-1, k+1, Lgrange_order), 1) =  j;
    
    dudy(map_ijk_to_idx(i,   j,   k,   Lgrange_order), 1) = -i + k;
    dudy(map_ijk_to_idx(i-1, j+1, k,   Lgrange_order), 1) = -i;
    dudy(map_ijk_to_idx(i-1, j,   k+1, Lgrange_order), 1) = -i;
    dudy(map_ijk_to_idx(i+1, j,   k-1, Lgrange_order), 1) =  k;
    dudy(map_ijk_to_idx(i,   j+1, k-1, Lgrange_order), 1) =  k;
    
    e = [dudx, dudy];
    end
    
    function idx = map_ijk_to_idx(i, j, k, n)
    idx = (n-i)*(n-i+1)/2 + (n-i-j) + 1;
    if i + j + k ~= n
        idx = -1;
    end
    if i<0 || j<0 || k<0
        idx = [];
    end
end


function ijk = create_ijk(n)
    ijk = zeros((n+1)*(n+2)/2, 3);
    index = 1;
    for p = n : -1 : 0
        for q = n-p : -1 : 0
            ijk(index, :) = [p, q, n-p-q];
            index = index + 1;
        end
    end
end

function [basis, nbasis] = Lagrange_basis(lagrange_order)
    n      = lagrange_order;
    nbasis = Lagrange_get_nbasis(lagrange_order);
    basis  = zeros(nbasis, 3);
    
    index = 1;
    %--------------- dof on verts ---------------
    basis(index, :) = [n, 0, 0];
    index = index + 1;
    if n > 0
        basis(index, :) = [0, n, 0];
        index = index + 1;
        basis(index, :) = [0, 0, n];
        index = index + 1;
    end
    %--------------- dof on edge 1 ---------------
    for p = n-1 : -1 : 1
        basis(index, :) = [0, p, n-p];
        index = index + 1;
    end
    %--------------- dof on edge 2 ---------------
    for p = n-1 : -1 : 1
        basis(index, :) = [n-p, 0, p];
        index = index + 1;
    end
    %--------------- dof on edge 3 ---------------
    for p = n-1 : -1 : 1
        basis(index, :) = [p, n-p, 0];
        index = index + 1;
    end
    %--------------- dof on element ---------------
    for p = n-2 : -1 : 1
        for q = n-1-p : -1 : 1
            basis(index, :) = [p, q, n-p-q];
            index = index + 1;
        end
    end
end

function nbasis = Lagrange_get_nbasis(lagrange_order)
    nbasis = (lagrange_order+1) * (lagrange_order+2) / 2;
end

function val = Lagrange_get_L1L2L3_ijk_value(i, j, k, x, y)
    val = (1-x-y)^i * x^j * y^k;
end


function y = Lagrange_integral_L1L2L3_ijk(i, j, k)
    y = I_intval(factorial(i) * factorial(j) * factorial(k)) / I_intval(factorial(i+j+k+2));
end

function y = Lagrange_integral_edge_L1L2L3_ijk(i, j, k, edge_idx)
    ijk = [i, j, k];
    if ijk(edge_idx) > 0
        y = 0;
    else
        y = I_intval(factorial(i) * factorial(j) * factorial(k)) / I_intval(factorial(i+j+k+1));
    end
end

function is_edge_bd = find_is_edge_bd(edge, bd_edge, ne, nb)
    edge_idx                = sort(edge,    2) * [ne; 1];
    bd_edge_idx             = sort(bd_edge, 2) * [ne; 1];
    [~, bd_edge_idx]        = ismember(bd_edge_idx, edge_idx);
    is_edge_bd              = zeros(ne, 1);
    is_edge_bd(bd_edge_idx) = ones(nb, 1);
end

function [tri, tri2edge] = find_tri2edge(edge, tri, ne, nt)
    tri2edge = zeros(nt, 3);
    edge_idx = sort(edge, 2) * [ne; 1];
    for k = 1:nt
        edge_local = [2 3; 1 3; 1 2];
        value = sort(reshape(tri(k, edge_local), 3, 2), 2) *[ne; 1];
        [~, idx] = ismember(value, edge_idx);
        tri2edge(k,:) = idx';
    end
end