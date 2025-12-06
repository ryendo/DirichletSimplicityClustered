<<<<<<< Updated upstream:VFEM2D_revised/Lagrange_dirichlet_eig_vectorized.m
function [eig_value, eig_func_no_bdry, eig_func_with_bdry, A_grad, A_L2, A_xx, A_xy, A_yy, bd_dof_idx] = Lagrange_dirichlet_eig_vectorized(Lagrange_order, vert, edge, tri, bd, neig)
=======
function [eig_value, eig_func_no_bdry, eig_func_with_bdry, A_grad, A_L2, A_xx, A_xy, A_yy, bd_dof_idx] = laplace_eig_lagrange_detailed(Lagrange_order, vert, edge, tri, bd, neig)
>>>>>>> Stashed changes:VFEM2D_revised/laplace_eig_lagrange_detailed.m

ne = size(edge, 1);
nt = size(tri,  1);
nb = size(bd,   1);

% if nt > 1000
%     disp('deal the mesh...');
% end
is_edge_bd = find_is_edge_bd(edge, bd, ne, nb);
[tri, tri2edge] = find_tri2edge(edge, tri, ne, nt);

% if nt > 1000
%     disp('create matrix...');
% end

[A_grad, A_L2, A_xx, A_xy, A_yy, bd_dof_idx, ndof] = Lagrange_dirichlet_eig_matrix(Lagrange_order, vert, edge, tri, tri2edge, is_edge_bd);

% if nt > 1000
%     disp('begin to solve eigenvalue problem');
% end
neig = min(neig, size(A_grad, 1));
[V, D] = eigs(I_mid(A_grad), I_mid(A_L2), neig, 'sm');

[eig_value, idx] = sort(diag(D));
eig_func_no_bdry = V(:, idx);
eig_func_with_bdry = I_intval(zeros(ndof,neig));
inside_dof_idx = 1:ndof; inside_dof_idx(bd_dof_idx)=[];
eig_func_with_bdry(inside_dof_idx,:) = eig_func_no_bdry;

dof_idx = 1:ndof;
dof_idx(bd_dof_idx)=[];
fem_func = I_intval(zeros(ndof,neig));
fem_func(dof_idx,:)=V(:, idx);

%eig_value = [(1:neig)', eig_value];
end

function [f_for_Hdiv_problem, point_for_Hdiv_problem] = Lagrange_get_femfunc_value_for_Hdiv_problem(Lagrange_order, vert, edge, tri, tri2edge, femfunc)
[basis, nbasis] = Lagrange_basis(Lagrange_order);
ne = size(edge, 1);
nt = size(tri, 1);
point_for_Hdiv_problem = cell(ne, 1);
point_for_Hdiv_problem_tri_idx = zeros(ne, 1);
edge_status = zeros(ne, 1);
for k = 1:nt
    edge_idx = tri2edge(k, :);
    for j = 1:3
        if edge_status(edge_idx(j)) > 0
            continue;
        else
            point = I_intval(zeros(Lagrange_order + 1, 2));
            point(1,   :) = vert(edge(edge_idx(j), 1), :);
            point(end, :) = vert(edge(edge_idx(j), 2), :);
            for i = 2:Lagrange_order
                point(i, :) = point(1, :) + (point(end, :)-point(1, :))*(i-1)/Lagrange_order;
            end
            point_for_Hdiv_problem{edge_idx(j), 1} = point;
            point_for_Hdiv_problem_tri_idx(edge_idx(j), 1) = k;
        end
    end
end

f_for_Hdiv_problem = I_intval(zeros(ne, Lagrange_order + 1));
for i = 1:ne
    point = point_for_Hdiv_problem{i, 1};
    tri_idx = point_for_Hdiv_problem_tri_idx(i);
    for j = 1:Lagrange_order+1        
        f_for_Hdiv_problem(i, j) = Lagrange_get_femfunc_value(Lagrange_order, basis, nbasis, ...
                                                              vert, edge, tri, tri2edge, ...
                                                              femfunc, point(j, 1), point(j, 2), tri_idx);
    end
end
end


function val = Lagrange_get_femfunc_value(Lagrange_order, basis, nbasis, vert, edge, tri, tri2edge, femfunc, x, y, tri_idx)
ne = size(edge, 1);
nv = size(vert, 1);
vert_idx = tri(tri_idx,:);
x1 = vert(vert_idx(1), 1); y1 = vert(vert_idx(1), 2);
x2 = vert(vert_idx(2), 1); y2 = vert(vert_idx(2), 2);
x3 = vert(vert_idx(3), 1); y3 = vert(vert_idx(3), 2);
B = [x2-x1, x3-x1; y2-y1, y3-y1];
Binv = [y3-y1, x1-x3; y1-y2, x2-x1] / det(B);
ksi_eta = Binv * [x-x1; y-y1];

edge_idx = tri2edge(tri_idx, :);
local_edge_start_vert = [2, 3, 1];
map_dof_idx_l2g = zeros(nbasis, 1);
map_dof_idx_l2g(1:3) = vert_idx;
for i = 1:3
    if edge(edge_idx(i), 1) == vert_idx(local_edge_start_vert(i))
        map_dof_idx_l2g(3+(Lagrange_order-1)*(i-1)+1:3+(Lagrange_order-1)*i) = nv+(Lagrange_order-1)*(edge_idx(i)-1)+1:1:nv+(Lagrange_order-1)*edge_idx(i);
    else
        map_dof_idx_l2g(3+(Lagrange_order-1)*(i-1)+1:3+(Lagrange_order-1)*i) = nv+(Lagrange_order-1)*edge_idx(i):-1:nv+(Lagrange_order-1)*(edge_idx(i)-1)+1;
    end
end
map_dof_idx_l2g(3+(Lagrange_order-1)*3+1:end) = nv + (Lagrange_order-1)*ne + ...
    ((Lagrange_order-1)*(Lagrange_order-2)/2*(tri_idx-1)+1:(Lagrange_order-1)*(Lagrange_order-2)/2*tri_idx);
coef = femfunc(map_dof_idx_l2g);
val = 0;
for i = 1:nbasis
    val = val + coef(i) * Lagrange_get_L1L2L3_ijk_value(basis(i, 1), basis(i, 2), basis(i, 3), ksi_eta(1), ksi_eta(2));
end
end


function [A_glob_grad, A_glob_L2, A_glob_ux_ux, A_glob_ux_uy, A_glob_uy_uy, bd_dof_idx, ndof] = Lagrange_dirichlet_eig_matrix(Lagrange_order, vert, edge, tri, tri2edge, is_edge_bd)
[basis, nbasis] = Lagrange_basis(Lagrange_order);
M_ip_elem  = Lagrange_inner_product_L1L2L3_all(Lagrange_order);
M_ip_edge1 = Lagrange_inner_product_edge_L1L2L3_all(Lagrange_order, 1);
M_ip_edge2 = Lagrange_inner_product_edge_L1L2L3_all(Lagrange_order, 2);
M_ip_edge3 = Lagrange_inner_product_edge_L1L2L3_all(Lagrange_order, 3);

M_ip_basis_grad_ijT_all = cell(nbasis, nbasis);
M_ip_basis_ij = zeros(nbasis, nbasis);
M_ip_basis_edge_ij_all = cell(3, 1);
M_ip_basis_edge1_ij = zeros(nbasis, nbasis);
M_ip_basis_edge2_ij = zeros(nbasis, nbasis);
M_ip_basis_edge3_ij = zeros(nbasis, nbasis);
for i = 1:nbasis
    for j = i:nbasis
        ei = Lagrange_create_coord_basis_grad(basis, i, Lagrange_order);
        ej = Lagrange_create_coord_basis_grad(basis, j, Lagrange_order);
        M_ip_basis_grad_ijT_all{j, i} = ei' * M_ip_elem * ej;
        
        ei = Lagrange_create_coord_basis(basis, i, Lagrange_order);
        ej = Lagrange_create_coord_basis(basis, j, Lagrange_order);
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

ndof = nv + (Lagrange_order-1)*ne + (Lagrange_order-1)*(Lagrange_order-2)/2*nt;

A_glob_grad = I_intval(zeros(ndof, ndof));
A_glob_L2 = I_intval(zeros(ndof, ndof));
A_glob_ux_ux = I_intval(zeros(ndof, ndof));
A_glob_ux_uy = I_intval(zeros(ndof, ndof));
A_glob_uy_uy = I_intval(zeros(ndof, ndof));
M = I_intval(zeros(ndof, ndof));

for k = 1:nt
    vert_idx = tri(k,:);
    x1 = vert(vert_idx(1), 1); y1 = vert(vert_idx(1), 2);
    x2 = vert(vert_idx(2), 1); y2 = vert(vert_idx(2), 2);
    x3 = vert(vert_idx(3), 1); y3 = vert(vert_idx(3), 2);
    B = [x2-x1, x3-x1; y2-y1, y3-y1];
    Binv = [y3-y1, x1-x3; y1-y2, x2-x1] / det(B);
    A_grad = I_intval(zeros(nbasis, nbasis));
    A_ux_vx = I_intval(zeros(nbasis, nbasis));
    A_ux_vy = I_intval(zeros(nbasis, nbasis));
    A_uy_vy = I_intval(zeros(nbasis, nbasis));
    for i = 1:nbasis
        for j = 1:i
            mat_base = Binv' * M_ip_basis_grad_ijT_all{i, j} * Binv;
            A_grad(i, j) = trace(mat_base) * det(B);
            A_ux_vx(i, j)   = mat_base(1,1) * det(B);
            A_ux_vy(i, j)   = (mat_base(1,2)+mat_base(2,1))/2 * det(B);
            A_uy_vy(i, j)   = mat_base(2,2) * det(B);
        end
    end
    A_grad = A_grad + tril(A_grad, -1)';
    A_ux_vx = A_ux_vx + tril(A_ux_vx, -1)';
    A_ux_vy = A_ux_vy + tril(A_ux_vy, -1)';
    A_uy_vy = A_uy_vy + tril(A_uy_vy, -1)';
    
    A_L2 = (M_ip_basis_ij + tril(M_ip_basis_ij, -1)') * det(B);
    
    edge_idx = tri2edge(k, :);
    edge_vec = vert(vert_idx([2,3,1]), :) - vert(vert_idx([3,1,2]), :);
    local_edge_start_vert = [2, 3, 1];
    map_dof_idx_l2g = zeros(nbasis, 1);
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
    
    A_glob_grad(map_dof_idx_l2g, map_dof_idx_l2g) = A_glob_grad(map_dof_idx_l2g, map_dof_idx_l2g) + A_grad;
    A_glob_ux_ux(map_dof_idx_l2g, map_dof_idx_l2g) = A_glob_ux_ux(map_dof_idx_l2g, map_dof_idx_l2g) + A_ux_vx;
    A_glob_ux_uy(map_dof_idx_l2g, map_dof_idx_l2g) = A_glob_ux_uy(map_dof_idx_l2g, map_dof_idx_l2g) + A_ux_vy;
    A_glob_uy_uy(map_dof_idx_l2g, map_dof_idx_l2g) = A_glob_uy_uy(map_dof_idx_l2g, map_dof_idx_l2g) + A_uy_vy;
    A_glob_L2(map_dof_idx_l2g, map_dof_idx_l2g) = A_glob_L2(map_dof_idx_l2g, map_dof_idx_l2g) + A_L2;

    for i = 1:3
        if is_edge_bd(edge_idx(i)) == 1
            M(map_dof_idx_l2g, map_dof_idx_l2g) = M(map_dof_idx_l2g, map_dof_idx_l2g) + sqrt(edge_vec(i, :)*edge_vec(i, :)') * M_ip_basis_edge_ij_all{i, 1};
        end
    end
end

bd_dof_idx=find(diag(M)>0);

A_glob_grad(bd_dof_idx,:)=[]; A_glob_grad(:,bd_dof_idx)=[];
A_glob_ux_ux(bd_dof_idx,:)=[]; A_glob_ux_ux(:,bd_dof_idx)=[];
A_glob_ux_uy(bd_dof_idx,:)=[]; A_glob_ux_uy(:,bd_dof_idx)=[];
A_glob_uy_uy(bd_dof_idx,:)=[]; A_glob_uy_uy(:,bd_dof_idx)=[];
A_glob_L2(bd_dof_idx,:)=[];   A_glob_L2(:,bd_dof_idx)=[];
end




function M_ip = Lagrange_inner_product_L1L2L3_all(Lagrange_order)
ijk = create_ijk(Lagrange_order);
len = size(ijk, 1);
M_ip = zeros(len, len);
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

function M_ip = Lagrange_inner_product_edge_L1L2L3_all(Lagrange_order, edge_idx)
ijk = create_ijk(Lagrange_order);
len = size(ijk, 1);
M_ip = zeros(len, len);
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

function val = Lagrange_get_L1L2L3_ijk_value(i, j, k, x, y)
val = (1-x-y)^i * x^j * y^k;
end


function y = Lagrange_integral_L1L2L3_ijk(i, j, k)
y = factorial(i) * factorial(j) * factorial(k) / factorial(i+j+k+2);
end

function y = Lagrange_integral_edge_L1L2L3_ijk(i, j, k, edge_idx)
ijk = [i, j, k];
if ijk(edge_idx) > 0
    y = 0;
else
    y = factorial(i) * factorial(j) * factorial(k) / factorial(i+j+k+1);
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