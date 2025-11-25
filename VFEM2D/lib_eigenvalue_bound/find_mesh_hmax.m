function hmax = find_mesh_hmax(vert, edges)
% FIND_MESH_HMAX  Compute the maximum edge length of a 2D mesh.
%
%   hmax = find_mesh_hmax(vert, edges)
%
% INPUT:
%   vert  - n×2 matrix. Each row is (x, y) coordinate of a vertex.
%   edges - m×2 matrix. Each row contains the vertex indices of an edge.
%
% OUTPUT:
%   hmax  - maximum edge length.

    % Extract coordinates of the first and second vertex of each edge
    v1 = vert(edges(:,1), :);   % m×2
    v2 = vert(edges(:,2), :);   % m×2

    % Compute squared distances for stability
    diff = v1 - v2;
    len_sq = diff(:,1).^2 + diff(:,2).^2;

    % Maximum edge length
    hmax = sqrt(I_intval(max(I_sup(len_sq))));
end