function eig_value = upper_eig_bound(Lagrange_order, vert, edge, tri, bd, neig)

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

[A_grad, A_L2] = Lagrange_dirichlet_eig_matrix_for_rho_vectorized(Lagrange_order, vert, edge, tri, tri2edge, is_edge_bd);

% if nt > 1000
%     disp(['begin to solve eigenvalue problem :size =',num2str(nt)]);
% end

[vec,approx_value] = eigs(I_mid(A_grad), I_mid(A_L2),neig,'sm');
LA = vec'*A_grad*vec;
LB = vec'*A_L2*vec;


global INTERVAL_MODE
if INTERVAL_MODE
  eig_value = veig(hull(LA,LA'), hull(LB,LB'), 1:neig);
else
  eig_value = eig(LA, LB);
end

eig_value = I_intval(I_sup(eig_value))';


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


function [A_glob_grad, A_glob_L2, A_glob_ux_ux, A_glob_ux_uy, A_glob_uy_uy, bd_dof_idx, ndof] = Lagrange_dirichlet_eig_matrix_for_rho_vectorized(Lagrange_order, vert, edge, tri, tri2edge, is_edge_bd)

% --- Precompute reference inner products ---
[basis, nbasis] = Lagrange_basis(Lagrange_order);
M_ip_elem  = Lagrange_inner_product_L1L2L3_all(Lagrange_order);
M_ip_edge1 = Lagrange_inner_product_edge_L1L2L3_all(Lagrange_order, 1);
M_ip_edge2 = Lagrange_inner_product_edge_L1L2L3_all(Lagrange_order, 2);
M_ip_edge3 = Lagrange_inner_product_edge_L1L2L3_all(Lagrange_order, 3);

% Build cell of local grad inner products and symmetric mass terms
M_ip_basis_grad = cell(nbasis, nbasis);
M_ip_basis_ij   = zeros(nbasis, nbasis);
M_ip_basis_edge = cell(3,1);
Mb1 = I_intval(zeros(nbasis, nbasis));
Mb2 = I_intval(zeros(nbasis, nbasis));
Mb3 = I_intval(zeros(nbasis, nbasis));
for i = 1:nbasis
    for j = i:nbasis
        eg = Lagrange_create_coord_basis_grad(basis, i, Lagrange_order);
        ej = Lagrange_create_coord_basis_grad(basis, j, Lagrange_order);
        M_ip_basis_grad{j,i} = eg' * M_ip_elem * ej;
        e  = Lagrange_create_coord_basis(basis, i, Lagrange_order);
        ej = Lagrange_create_coord_basis(basis, j, Lagrange_order);
        M_ip_basis_ij(j,i)   = e' * M_ip_elem * ej;
        Mb1(j,i) = e' * M_ip_edge1 * ej;
        Mb2(j,i) = e' * M_ip_edge2 * ej;
        Mb3(j,i) = e' * M_ip_edge3 * ej;
    end
end
for i = 1:nbasis
    for j = i+1:nbasis
        M_ip_basis_grad{i,j} = M_ip_basis_grad{j,i}';
    end
end

Mb1 = Mb1 + tril(Mb1, -1)';
Mb2 = Mb2 + tril(Mb2, -1)';
Mb3 = Mb3 + tril(Mb3, -1)';

M_ip_basis_edge = { Mb1, Mb2, Mb3 };


% --- Degrees of freedom count ---
nv = size(vert,1);
ne = size(edge,1);
nt = size(tri,1);
ndof = nv + (Lagrange_order-1)*ne + (Lagrange_order-1)*(Lagrange_order-2)/2*nt;

% --- Build gradient S_vec (4×nbasis^2) ---
S_vec = I_zeros(4, nbasis^2);
cnt = 1;
for j = 1:nbasis
    for i = 1:nbasis
        Mij = M_ip_basis_grad{i,j};
        S_vec(:,cnt) = [Mij(1,1); Mij(2,1); Mij(1,2); Mij(2,2)];
        cnt = cnt + 1;
    end
end

% --- Precompute geometry transforms in batch ---
X1 = vert(tri(:,1),:);
X2 = vert(tri(:,2),:);
X3 = vert(tri(:,3),:);

% B(:,:,k) and its det & inverse
B    = I_intval(zeros(2,2,nt));
B(1,1,:) = X2(:,1) - X1(:,1);
B(2,1,:) = X2(:,2) - X1(:,2);
B(1,2,:) = X3(:,1) - X1(:,1);
B(2,2,:) = X3(:,2) - X1(:,2);

detB = squeeze(B(1,1,:).*B(2,2,:) - B(1,2,:).*B(2,1,:));
Binv = I_intval(zeros(2,2,nt));
Binv(1,1,:) =  B(2,2,:);
Binv(1,2,:) = -B(1,2,:);
Binv(2,1,:) = -B(2,1,:);
Binv(2,2,:) =  B(1,1,:);
for k = 1:nt
    Binv(:,:,k) = Binv(:,:,k) / detB(k);
end

% Precompute K(:,:,k) = kron(Binv(:,:,k)', Binv(:,:,k)')
K = I_intval(zeros(4,4,nt));
for k = 1:nt
    K(:,:,k) = kron(Binv(:,:,k)', Binv(:,:,k)');
end

% Compute M_base_all (4×nb^2×nt)
p = size(S_vec,2);
M_tmp = I_intval(zeros(4,p,nt));
for m = 1:4
    Km = reshape(K(:,m,:),4,1,nt);
    Sm = reshape(S_vec(m,:),1,p,1);
    M_tmp = M_tmp + Km .* Sm;
end
M_base_all = M_tmp .* reshape(detB,1,1,nt);

% --- Initialize global (intval) matrices ---
mat0 = zeros(ndof, ndof);
mat_ival = I_intval(mat0);
A_glob_grad  = mat_ival;
A_glob_L2    = mat_ival;
M_bd        = mat_ival;

% --- Assemble via element loop ---
for k = 1:nt
    % local-to-global DOF mapping
    vert_idx = tri(k,:);
    map = zeros(nbasis,1);
    map(1:3) = vert_idx;
    % edge DOFs
    for e_id = 1:3
        e = tri2edge(k,e_id);
        v0 = vert_idx(mod(e_id,3)+1);
        r  = 3+(Lagrange_order-1)*(e_id-1)+1 : 3+(Lagrange_order-1)*e_id;
        if edge(e,1)==v0
            map(r) = nv + (Lagrange_order-1)*(e-1) + (1:(Lagrange_order-1));
        else
            map(r) = nv + (Lagrange_order-1)*e - (0:(Lagrange_order-1)-1);
        end
    end
    % interior DOFs
    base = nv + (Lagrange_order-1)*ne;
    n_int = (Lagrange_order-1)*(Lagrange_order-2)/2;
    map(3+(Lagrange_order-1)*3+1:end) = base + (k-1)*n_int + (1:n_int);

    % extract precomputed M_base_vec
    Mbv = M_base_all(:,:,k);
    A_trace = Mbv(1,:) + Mbv(4,:);

    % reshape & symmetrize
    Ag = reshape(A_trace, nbasis, nbasis);  Ag  = tril(Ag)  + tril(Ag,-1)';
    AL2 = (M_ip_basis_ij + tril(M_ip_basis_ij,-1)') * detB(k);

    % add local to global
    A_glob_grad(map,map)  = A_glob_grad(map,map)  + Ag;
    A_glob_L2(map,map)    = A_glob_L2(map,map)    + AL2;

    % boundary edge contributions
    eids = tri2edge(k,:);
    evs  = vert(vert_idx([2,3,1]),:) - vert(vert_idx([3,1,2]),:);
    for eid = 1:3
        if is_edge_bd(eids(eid))
            len = sqrt(evs(eid,:)*evs(eid,:)');
            M_bd(map,map) = M_bd(map,map) + len * M_ip_basis_edge{eid};
        end
    end
end

% --- Extract Dirichlet DOFs & eliminate ---
bd_dof_idx = find(diag(M_bd)>0);
A_glob_grad(bd_dof_idx,:) = [];
A_glob_grad(:,bd_dof_idx) = [];
A_glob_L2(bd_dof_idx,:)   = [];
A_glob_L2(:,bd_dof_idx)   = [];
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
