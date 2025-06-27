function [x, FA, F, Fstiff] = RT_Hdiv_problem(RT_order, vert, edge, tri, bd, f)

% dim of f: ne * (RT_order+1)
ne = size(edge, 1);
nt = size(tri,  1);
nb = size(bd,   1);

if nt > 1000
    disp('deal the mesh...');
end
is_edge_bd = find_is_edge_bd(edge, bd, ne, nb);
[tri, tri2edge] = find_tri2edge(edge, tri, ne, nt);

if nt > 1000
    disp('create matrix...');
end
A = RT_Hdiv_stiff_matrix(RT_order, vert, edge, tri, tri2edge);
F = zeros(size(A, 1), 1);
Fstiff = full(A);

if nt > 1000
    disp('deal with boundary condition...');
end
[A, F] = RT_Hdiv_boundary_condition(A, F, RT_order, vert, edge, tri, tri2edge, is_edge_bd, f);
FA = full(A);

if nt > 1000
    disp('solve fem system...');
end
x = A\F;
end

function [A, F] = RT_Hdiv_boundary_condition(A, F, RT_order, vert, edge, tri, tri2edge, is_edge_bd, f)
f_coord = f;
if RT_order == 1
    f_coord = f;
elseif RT_order == 2
    f_coord(:, 2) = -(f(:, 1) - 4*f(:, 2) + f(:, 3));
else
    error('RT_order can only be 1 or 2!');
end

nt = size(tri, 1);
for k = 1:nt
    vert_idx = tri(k,:);
    edge_vec = vert(vert_idx([2,3,1]), :) - vert(vert_idx([3,1,2]), :);
    edge_idx = tri2edge(k, :);
    local_edge_start_vert = [2, 3, 1];
    for i = 1:3
        if is_edge_bd(edge_idx(i)) == 1
            if edge(edge_idx(i), 1) == vert_idx(local_edge_start_vert(i))
                map_dof_idx_l2g = (RT_order+1)*(edge_idx(i)-1)+1:1:(RT_order+1)*edge_idx(i);
                F(map_dof_idx_l2g) = f_coord(edge_idx(i), 1:RT_order+1);
            else
                map_dof_idx_l2g = (RT_order+1)*edge_idx(i):-1:(RT_order+1)*(edge_idx(i)-1)+1;
                F(map_dof_idx_l2g) = -f_coord(edge_idx(i), RT_order+1:-1:1);
            end
            
            A(map_dof_idx_l2g, :) = 0;
            A(map_dof_idx_l2g, map_dof_idx_l2g) =  eye(RT_order+1) / sqrt(edge_vec(i, :)*edge_vec(i, :)');
        end
    end
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

function A = RT_Hdiv_stiff_matrix(RT_order, vert, edge, tri, tri2edge)
[basis_abc, basis_ijk, nbasis] = RT_basis(RT_order);
Mnp1 = RT_inner_product_L1L2L3_all(RT_order + 1);
Mn   = RT_inner_product_L1L2L3_all(RT_order);

M_ip_basis_ijT_all = cell(nbasis);
M_ip_divbasis = zeros(nbasis, nbasis);
for i = 1:nbasis
    for j = i:nbasis
        ei = RT_create_coord_basis(basis_abc, basis_ijk, i, RT_order);
        ej = RT_create_coord_basis(basis_abc, basis_ijk, j, RT_order);
        M_ip_basis_ijT_all{j, i} = ei' * Mnp1 * ej;
        
        ei = RT_create_coord_basis_div(basis_abc, basis_ijk, i, RT_order);
        ej = RT_create_coord_basis_div(basis_abc, basis_ijk, j, RT_order);
        M_ip_divbasis(j, i) = ei' * Mn * ej;
    end
end

ne = size(edge, 1);
nt = size(tri,  1);

ndof = (RT_order+1)*ne + RT_order*(RT_order+1)*nt;
A = sparse(ndof, ndof);

for k = 1:nt
    vert_idx = tri(k,:);
    x1 = vert(vert_idx(1), 1); y1 = vert(vert_idx(1), 2);
    x2 = vert(vert_idx(2), 1); y2 = vert(vert_idx(2), 2);
    x3 = vert(vert_idx(3), 1); y3 = vert(vert_idx(3), 2);
    B = [x2-x1, x3-x1; y2-y1, y3-y1];
    A1 = zeros(nbasis, nbasis);
    for i = 1:nbasis
        for j = 1:i
            A1(i, j) = trace(B * M_ip_basis_ijT_all{i, j} * B') / det(B);
        end
    end
    A1 = A1 + tril(A1, -1)';
    
    A2 = (M_ip_divbasis + tril(M_ip_divbasis, -1)') / det(B);
    
    edge_idx = tri2edge(k, :);
    local_edge_start_vert = [2, 3, 1];
    map_dof_idx_l2g = zeros(nbasis, 1);
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
    A(map_dof_idx_l2g, map_dof_idx_l2g) = A(map_dof_idx_l2g, map_dof_idx_l2g) + diag(P)*(A1+A2)*diag(P);
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
nbasis = (RT_order+1) * (RT_order+3);
end

function y = RT_integral_L1L2L3_ijk(i, j, k)
y = factorial(i) * factorial(j) * factorial(k) / factorial(i+j+k+2);
end