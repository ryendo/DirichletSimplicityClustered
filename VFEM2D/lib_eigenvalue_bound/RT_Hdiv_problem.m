function [mat_b_w_w] = RT_Hdiv_problem(mesh, RT_order, f)
    % dim of f: ne * (RT_order+1)
    global INTERVAL_MODE

    vert = I_intval(mesh.nodes);
    edge = mesh.edges;
    tri  = mesh.elements;

    nt = size(tri,  1);
    ne = size(edge, 1);
    [tri, tri2edge] = find_tri2edge(edge, tri, ne, nt);
    
    if nt > 1000
        disp('create matrix...');
    end
    % A = RT_Hdiv_stiff_matrix(RT_order, vert, edge, tri, tri2edge);
    [A,B] = RT_mixed_matrix(RT_order, vert, edge, tri,tri2edge);
    disp("Size of matrix for RT:")
    disp(size(A,1))    
    
    if nt > 1000
        disp('deal with boundary condition...');
    end
    F = RT_mixed_set_f(RT_order, vert, edge, tri, tri2edge, f);

    % [A, F] = RT_Hdiv_boundary_condition(A, F, RT_order, vert, edge, tri, tri2edge, is_edge_bd, f);
    if INTERVAL_MODE
        %% Improve computation for mat_b_w_w
        if nt > 1000
            disp('solve RT fem system...');
        end
        H = B'*(I_solve(A,B));
        Y = I_solve(H,-F);
        mat_b_w_w = Y'*(-F);

    else
    %% Stardard way to calculate mat_b_w_w (w: Goerisch w term)
        m = size(B,2);
        n = size(A,1);
        mat_mixed = sparse([A,B;B', zeros(m,m)]);
        f_mixed = I_zeros(size(mat_mixed,1), size(f,2));
        f_mixed(n+1:end,:) = - F;
        if nt > 1000
            disp('solve RT fem system...');
        end
        tic;
        x = I_solve(mat_mixed,f_mixed);
        toc;
        RT_vec = x(1:n,:);
        mat_b_w_w = RT_vec' * A * RT_vec;

    end

end

function F = RT_mixed_set_f(RT_order, vert, edge, tri, tri2edge, f)
    Lagrange_order = RT_order;

    dof_dg_elt = RT_get_Bernstein_polynomial_nbasis(RT_order);
    [basis, n_dg_elt] = Lagrange_basis(Lagrange_order);
    Mn   = RT_inner_product_L1L2L3_all(RT_order);

    ne = size(edge, 1);
    nt = size(tri,  1);
    nv = size(vert,1);

    M_ip_L2 = zeros(n_dg_elt, n_dg_elt);
    for i = 1:n_dg_elt
        ei = Lagrange_create_coord_basis(basis, i, Lagrange_order);
        for j = 1:n_dg_elt            
            ej = Lagrange_create_coord_basis(basis, j, Lagrange_order);
            M_ip_L2(i, j) = ei' * Mn * ej;
        end
    end
    
    
    n_dim_f = size(f,2);
    n_dg = nt*n_dg_elt; 
    F = I_zeros(n_dg,n_dim_f);
    for k = 1:nt
        vert_idx = tri(k,:);
        edge_idx = tri2edge(k, :);
        x1 = vert(vert_idx(1), 1); y1 = vert(vert_idx(1), 2);
        x2 = vert(vert_idx(2), 1); y2 = vert(vert_idx(2), 2);
        x3 = vert(vert_idx(3), 1); y3 = vert(vert_idx(3), 2);
        B = [x2-x1, x3-x1; y2-y1, y3-y1];
        
        %[Liu] edge 1: two vertices as [3,2] or [2,3]
        %[Liu] edge 2: two vertices as [1,3] or [3,1]
        %[Liu] edge 1: two vertices as [2,1] or [1,2]
        local_edge_start_vert = [2, 3, 1];

        map_dof_idx_l2g_dg = (k-1)*dof_dg_elt+(1:dof_dg_elt);


        map_dof_idx_l2g = zeros(dof_dg_elt, 1);
        map_dof_idx_l2g(1:3) = vert_idx;
        for i = 1:3
            if edge(edge_idx(i), 1) == vert_idx(local_edge_start_vert(i))
                map_dof_idx_l2g(3+(Lagrange_order-1)*(i-1)+1:3+(Lagrange_order-1)*i) = nv+(Lagrange_order-1)*(edge_idx(i)-1)+1:1:nv+(Lagrange_order-1)*edge_idx(i);
            else
                map_dof_idx_l2g(3+(Lagrange_order-1)*(i-1)+1:3+(Lagrange_order-1)*i) = nv+(Lagrange_order-1)*edge_idx(i):-1:nv+(Lagrange_order-1)*(edge_idx(i)-1)+1;
            end
            
        end
        map_dof_idx_l2g(3+(Lagrange_order-1)*3+1:end) = nv + (Lagrange_order-1)*ne + ...
                                                        ((Lagrange_order-1)*(Lagrange_order-2)/2*(k-1)+1:(Lagrange_order-1)*(Lagrange_order-2)/2*k);
        
        F(map_dof_idx_l2g_dg,:) = M_ip_L2*f(map_dof_idx_l2g,:) * det(B);
    end
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


function [mat_A,mat_B] = RT_mixed_matrix(RT_order, vert, edge, tri,tri2edge)
    Lagrange_order = RT_order;
    [basis, n_dg_elt] = Lagrange_basis(Lagrange_order);
    [basis_abc, basis_ijk, nbasis] = RT_basis(RT_order);
    Mnp1 = RT_inner_product_L1L2L3_all(RT_order + 1);%[order to plus one]
    Mn   = RT_inner_product_L1L2L3_all(RT_order);


    ne = size(edge, 1);
    nt = size(tri,  1);
    n_dg = nt*n_dg_elt; 

    M_ip_basis_ijT_all = cell(nbasis);
    M_ip_div_dg = zeros(nbasis, n_dg_elt);
    for i = 1:nbasis
        ei = RT_create_coord_basis(basis_abc, basis_ijk, i, RT_order);
        for j = i:nbasis
            ej = RT_create_coord_basis(basis_abc, basis_ijk, j, RT_order);
            M_ip_basis_ijT_all{j, i} = ei' * Mnp1 * ej;
        end
    end
    for i = 1:nbasis
        ei = RT_create_coord_basis_div(basis_abc, basis_ijk, i, RT_order);
        for j = 1:n_dg_elt            
            ej = Lagrange_create_coord_basis(basis, j, Lagrange_order);
            M_ip_div_dg(i, j) = ei' * Mn * ej;
        end
    end
    
    
    %[Liu] dof of RT elements. In case RT_order>1, there are interior DOFs
    ndof = (RT_order+1)*ne + RT_order*(RT_order+1)*nt;
    mat_A = I_sparse(ndof, ndof);
    mat_B = I_sparse(ndof, n_dg);
    
    for k = 1:nt
        vert_idx = tri(k,:);
        x1 = vert(vert_idx(1), 1); y1 = vert(vert_idx(1), 2);
        x2 = vert(vert_idx(2), 1); y2 = vert(vert_idx(2), 2);
        x3 = vert(vert_idx(3), 1); y3 = vert(vert_idx(3), 2);
        B = [x2-x1, x3-x1; y2-y1, y3-y1];
        A1 = I_zeros(nbasis, nbasis);
        for i = 1:nbasis
            for j = 1:i
                A1(i, j) = trace(B * M_ip_basis_ijT_all{i, j} * B') / det(B);
            end
        end
        A1 = A1 + tril(A1, -1)';
        
        A2 = M_ip_div_dg ; % [Liu] no scaling of det(B) is needed here. sqrt(det(B)) ;% /det(B) 
        
        edge_idx = tri2edge(k, :);
        local_edge_start_vert = [2, 3, 1];
        map_dof_idx_l2g = zeros(nbasis, 1);
        map_dof_idx_l2g_dg = (k-1)*n_dg_elt + (1:n_dg_elt);
        P = ones(nbasis, 1);
        for i = 1:3
            if edge(edge_idx(i), 1) == vert_idx(local_edge_start_vert(i))
                map_dof_idx_l2g((RT_order+1)*(i-1)+1:(RT_order+1)*i) = (RT_order+1)*(edge_idx(i)-1)+1:1:(RT_order+1)*edge_idx(i);
            else
                map_dof_idx_l2g((RT_order+1)*(i-1)+1:(RT_order+1)*i) = (RT_order+1)*edge_idx(i):-1:(RT_order+1)*(edge_idx(i)-1)+1;
                P((RT_order+1)*(i-1)+1:(RT_order+1)*i) = -1;
            end
        end
        map_dof_idx_l2g((RT_order+1)*3+1:end) = (RT_order+1)*ne + (RT_order*(RT_order+1)*(k-1)+1:RT_order*(RT_order+1)*k);
        mat_A(map_dof_idx_l2g, map_dof_idx_l2g) = mat_A(map_dof_idx_l2g, map_dof_idx_l2g) + diag(P)*(A1)*diag(P);
        mat_B(map_dof_idx_l2g, map_dof_idx_l2g_dg) = diag(P)*A2;
    end
end



function M = RT_inner_product_L1L2L3_all(RT_order)
    ijk = create_ijk(RT_order);
    len = size(ijk, 1);
    M = zeros(len, len);
    for p = 1:len
        for q = p:len
            pi = ijk(p, 1);
            pj = ijk(p, 2);
            pk = ijk(p, 3);
            qi = ijk(q, 1);
            qj = ijk(q, 2);
            qk = ijk(q, 3);
            M(p, q) = RT_integral_L1L2L3_ijk(pi+qi, pj+qj, pk+qk);
        end
    end
    M = M + triu(M, 1)';
end

function e = RT_create_coord_basis(basis_abc, basis_ijk, idx, RT_order)
    len = RT_get_Bernstein_polynomial_nbasis(RT_order + 1);
    e = zeros(len, 2);
    
    a = basis_abc(idx, 1);
    b = basis_abc(idx, 2);
    c = basis_abc(idx, 3);
    i = basis_ijk(idx, 1);
    j = basis_ijk(idx, 2);
    k = basis_ijk(idx, 3);
    
    e(RT_map_ijk_to_idx(i+1, j,   k,   RT_order+1), 1) = a;
    e(RT_map_ijk_to_idx(i,   j+1, k,   RT_order+1), 1) = a+c;
    e(RT_map_ijk_to_idx(i,   j,   k+1, RT_order+1), 1) = a;
    
    e(RT_map_ijk_to_idx(i+1, j,   k,   RT_order+1), 2) = b;
    e(RT_map_ijk_to_idx(i,   j+1, k,   RT_order+1), 2) = b;
    e(RT_map_ijk_to_idx(i,   j,   k+1, RT_order+1), 2) = b+c;
end

function e = RT_create_coord_basis_div(basis_abc, basis_ijk, idx, RT_order)
len = RT_get_Bernstein_polynomial_nbasis(RT_order);
e = zeros(len, 1);
f = zeros(len, 1);

a = basis_abc(idx, 1);
b = basis_abc(idx, 2);
c = basis_abc(idx, 3);
i = basis_ijk(idx, 1);
j = basis_ijk(idx, 2);
k = basis_ijk(idx, 3);

e(RT_map_ijk_to_idx(i,   j,   k,   RT_order)) = a*(j-i) + c*(j+1);
e(RT_map_ijk_to_idx(i-1, j+1, k,   RT_order)) = -(a+c) * i;
e(RT_map_ijk_to_idx(i-1, j,   k+1, RT_order)) = -a * i;
e(RT_map_ijk_to_idx(i+1, j-1, k,   RT_order)) = a * j;
e(RT_map_ijk_to_idx(i,   j-1, k+1, RT_order)) = a * j;

f(RT_map_ijk_to_idx(i,   j,   k,   RT_order)) = b*(k-i) + c*(k+1);
f(RT_map_ijk_to_idx(i-1, j,   k+1, RT_order)) = -(b+c) * i;
f(RT_map_ijk_to_idx(i-1, j+1, k,   RT_order)) = -b * i;
f(RT_map_ijk_to_idx(i+1, j  , k-1, RT_order)) = b * k;
f(RT_map_ijk_to_idx(i,   j+1, k-1, RT_order)) = b * k;

e = e+f;
end

function idx = RT_map_ijk_to_idx(i, j, k, n)
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

function [basis_abc, basis_ijk, nbasis] = RT_basis(RT_order)
%RT_order = 2;
n = RT_order;
nbasis = RT_get_nbasis(RT_order);
basis_abc = zeros(nbasis, 3);
basis_ijk = zeros(nbasis, 3);

index = 1;

%--------------- edge 1 ---------------
basis_abc(index, :) = [0, 0, 1];
basis_ijk(index, :) = [0, n, 0];
index = index + 1;
for p = n-1 : -1 : 1
    basis_abc(index, :) = [1, 0, 0];
    basis_ijk(index, :) = [0, p, n-p];
    index = index + 1;
end
if n~= 0
    basis_abc(index, :) = [0, 1, 0];
    basis_ijk(index, :) = [0, 0, n];
    index = index + 1;
end

%--------------- edge 2 ---------------
basis_abc(index, :) = [-1, 0, 1];
basis_ijk(index, :) = [0, 0, n];
index = index + 1;
for p = n-1 : -1 : 1
    basis_abc(index, :) = [-1, 0, 0];
    basis_ijk(index, :) = [n-p, 0, p];
    index = index + 1;
end
if n~= 0
    basis_abc(index, :) = [-1, 0, 0];
    basis_ijk(index, :) = [n, 0, 0];
    index = index + 1;
end

%--------------- edge 3 ---------------
basis_abc(index, :) = [0, -1, 1];
basis_ijk(index, :) = [n, 0, 0];
index = index + 1;
for p = n-1 : -1 : 1
    basis_abc(index, :) = [0, -1, 0];
    basis_ijk(index, :) = [p, n-p, 0];
    index = index + 1;
end
if n~= 0
    basis_abc(index, :) = [0, -1, 1];
    basis_ijk(index, :) = [0, n, 0];
    index = index + 1;
end


if n >= 1
    basis_abc(index, :) = [0, 0, 1];
    basis_ijk(index, :) = [n, 0, 0];
    index = index + 1;
    basis_abc(index, :) = [-1, 0, 1];
    basis_ijk(index, :) = [0, n, 0];
    index = index + 1;
end

if n >= 2
    for p = n-1 : -1 : 1
        basis_abc(index, :) = [1, -1, 0];
        basis_ijk(index, :) = [0, p, n-p];
        index = index + 1;
    end
    
    for p = n-1 : -1 : 1
        basis_abc(index, :) = [0, 1, 0];
        basis_ijk(index, :) = [n-p, 0, p];
        index = index + 1;
    end
    
    for p = n-1 : -1 : 1
        basis_abc(index, :) = [-1, 0, 0];
        basis_ijk(index, :) = [p, n-p, 0];
        index = index + 1;
    end
    
    for p = n-1 : -1 : 1
        basis_abc(index, :) = [0, 0, 1];
        basis_ijk(index, :) = [p, n-p, 0];
        index = index + 1;
    end
end

if n >= 3
    for p = n-2 : -1 : 1
        for q = n-1-p : -1 : 1
            basis_abc(index, :) = [1, 0, 0];
            basis_ijk(index, :) = [p, q, n-p-q];
            index = index + 1;
        end
    end
    
    for p = n-2 : -1 : 1
        for q = n-1-p : -1 : 1
            basis_abc(index, :) = [0, 1, 0];
            basis_ijk(index, :) = [p, q, n-p-q];
            index = index + 1;
        end
    end
end

end

function nbasis = RT_get_Bernstein_polynomial_nbasis(degree)
    nbasis = (degree+1) * (degree+2) / 2;
end

function nbasis = RT_get_nbasis(RT_order)
    %[Liu] Minimal order is zero. 
    nbasis = (RT_order+1) * (RT_order+3);
end

function y = RT_integral_L1L2L3_ijk(i, j, k)
y = factorial(i) * factorial(j) * factorial(k) / factorial(i+j+k+2);
end

%%%%%%% For lagrange element %%%%%%%%

function e = Lagrange_create_coord_basis(basis, idx, Lgrange_order)
len = Lagrange_get_nbasis(Lgrange_order);
e = zeros(len, 1);

i = basis(idx, 1);
j = basis(idx, 2);
k = basis(idx, 3);
e(map_ijk_to_idx(i, j, k, Lgrange_order), 1) = 1;
end




function ijk = LG_create_ijk(n)
ijk = zeros((n+1)*(n+2)/2, 3);
index = 1;
for p = n : -1 : 0
    for q = n-p : -1 : 0
        ijk(index, :) = [p, q, n-p-q];
        index = index + 1;
    end
end
end

function [basis, nbasis] = Lagrange_basis(Lagrange_order)
n      = Lagrange_order;
nbasis = Lagrange_get_nbasis(Lagrange_order);
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

function nbasis = Lagrange_get_nbasis(Lagrange_order)
    nbasis = (Lagrange_order+1) * (Lagrange_order+2) / 2;
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