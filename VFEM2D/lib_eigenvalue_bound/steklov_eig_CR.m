function [eig_value, Ch] = steklov_eig_CR(vert, edge, tri, bd_edge, neig)
%
% Calculate the eigenvalues of Steklov eigenvalue problem
% using Crouzeix-Raviart element
%
%                  -\Delta u + u = 0            in \Omega
% \frac{\partial u}{\partial\nu} = \lambda u    on \partial \Omega
% where \nu denotes the outer normal direction
%
% The variational problem is
% \int_{\Omega} (\nabla u \nabla v + uv) d\Omega = \lambda\int_{\partial\Omega} uv ds
%

ne = size(edge,    1);
nt = size(tri,     1);
nb = size(bd_edge, 1);

if nt > 1000
    disp('deal the mesh...');
end
is_edge_bd = find_is_edge_bd(edge, bd_edge, ne, nb);
[tri, tri2edge] = find_tri2edge(edge, tri, ne, nt);

if nt > 1000
    disp('create matrix...');
end
[A, M] = create_matrix(vert, edge, tri, tri2edge, is_edge_bd);

if nt > 1000
    disp('begin to solve eigenvalue problem');
end
[~, d] = eigs(A, M, min(neig, ne), 'sm');
eig_value = sort(diag(d));

if nt > 1000
    disp('begin to compute Ch');
end
Ch = compute_Ch(vert, tri, tri2edge, is_edge_bd, eig_value);
end


%=========================================================================%
%============================  SUB FUNCTIONS  ============================%
%=========================================================================%
function is_edge_bd = find_is_edge_bd(edge, bd_edge, ne, nb)
edge_idx    = sort(edge,    2) * [ne; 1];
bd_edge_idx = sort(bd_edge, 2) * [ne; 1];

[~, bd_edge_idx] = ismember(bd_edge_idx, edge_idx);

is_edge_bd  = zeros(ne, 1);
is_edge_bd(bd_edge_idx) = ones(nb, 1);
end

%-------------------------------------------------------------------------%
function [tri, tri2edge] = find_tri2edge(edge, tri, ne, nt)
tri2edge = zeros(nt, 3);
edge_idx = sort(edge, 2) * [ne; 1];
for k = 1:nt
    edge_local = [2 3; 1 3; 1 2];
    value = sort(reshape(tri(k, edge_local), 3, 2), 2) *[ne; 1];
    [~, idx] = ismember(value, edge_idx);
    tri2edge(k,:) = idx';
    %
    %     %promise the 3rd edge is boundary edge
    %     if is_edge_bd(tri2edge(k, 1)) == 1
    %         tri(k, :) = [tri(k, 2), tri(k, 3), tri(k, 1)];
    %         tri2edge(k,:) = [tri2edge(k, 2), tri2edge(k, 3), tri2edge(k, 1)];
    %     elseif is_edge_bd(tri2edge(k, 2)) == 1
    %         tri(k, :) = [tri(k, 3), tri(k, 1), tri(k, 2)];
    %         tri2edge(k,:) = [tri2edge(k, 3), tri2edge(k, 1), tri2edge(k, 2)];
    %     end
end
end

%-------------------------------------------------------------------------%
function [A, M] = create_matrix(vert, edge, tri, tri2edge, is_edge_bd)
nt = size(tri, 1);
ne = size(edge,1);

A  = sparse(ne, ne);
M  = sparse(ne, ne);

for k=1:nt
    edge_idx = tri2edge(k, :);
    vert_idx = tri(k,:);
    
    edge_vec = vert(vert_idx([2,3,1]), :) - vert(vert_idx([3,1,2]), :);
    S = 0.5 * det(edge_vec(2:3, :));
    e_e = edge_vec * edge_vec';
    
    A(edge_idx, edge_idx) = A(edge_idx, edge_idx) + S*eye(3,3)/3.0 + e_e/S;
    
    if is_edge_bd(edge_idx(1)) == 1
        int_b = [1, 0, 0; 0, 1/3, -1/3; 0, -1/3, 1/3];
        M(edge_idx, edge_idx) = M(edge_idx, edge_idx) + sqrt(edge_vec(1, :)*edge_vec(1, :)') * int_b;
    end
    if is_edge_bd(edge_idx(2)) == 1
        int_b = [1/3, 0, -1/3; 0, 1, 0; -1/3, 0, 1/3];
        M(edge_idx, edge_idx) = M(edge_idx, edge_idx) + sqrt(edge_vec(2, :)*edge_vec(2, :)') * int_b;
    end
    if is_edge_bd(edge_idx(3)) == 1
        int_b = [1/3, -1/3, 0; -1/3, 1/3, 0; 0, 0, 1;];
        M(edge_idx, edge_idx) = M(edge_idx, edge_idx) + sqrt(edge_vec(3, :)*edge_vec(3, :)') * int_b;
    end
end
end

%-------------------------------------------------------------------------%
function Ch = compute_Ch(vert, tri, tri2edge, is_edge_bd, eig_value)
C1 = 0;
C2 = 0;
for k = 1:size(tri,1)
    edge_idx = tri2edge(k,:);
    vert_idx = tri(k,:);
    edge_vec = vert(vert_idx([2,3,1]), :) - vert(vert_idx([3,1,2]), :);
    edge_len = sqrt(edge_vec(:,1).^2 + edge_vec(:,2).^2);
    
    HK = max(edge_len);
    C2 = max(C2, HK);
    
    if sum(is_edge_bd(edge_idx)) > 0
        [hK, max_ind] = max(edge_len(1:2));
        cos_theta = (edge_vec(max_ind) * edge_vec(3)') / (edge_len(max_ind) * edge_len(3));
        C1 = max(C1, HK / sqrt(hK * sqrt(1 - cos_theta^2)));
    end
end
Ch = (1 + sqrt(2) * 0.1893) * C1 + 0.1893 * C2 / sqrt(eig_value(1));

Ch = [Ch, C1, C2];
end