function mesh = read_mesh_from_folder(folder_name)
%% Load data from specified folder to create a structure mesh.
%   mesh = read_mesh_from_folder(folder_name)
%
% Inputs:
%   folder_name - The folder for mesh data.
%
% Outputs:
%   mesh - Structure (fields: nodes, edges, elements, domain, min_edge_length)
%          nodes: Each row contains [x, y] coordinates (z component in XML is ignored)
%          elements: Each row contains 3 vertex indices (triangle elements in XML; Note: XML is 0-based, so +1 is added for MATLAB)
%          edges: Unique edges extracted from elements (each edge is sorted in ascending order)
%          domain: Stores vertices of the first triangle in [1, 3, 2] order as a sample
%          bd_edges: Boundary edges represented by two end node indices.
%          %min_edge_length: The minimum value among the maximum edge lengths of each triangle element
%
% Example:
%   mesh = read_dolfin_mesh('./mesh8/');
% Revision history:
%   2024/06/10: Xuefeng Liu, First edition.

    mesh_path = folder_name;
    vert = load([mesh_path, 'vert.dat']);
    edges = load([mesh_path, 'edge.dat']);
    tri  = load([mesh_path, 'tri.dat']);
    bd_edges   = load([mesh_path, 'bd.dat']);
    
    domain = [];

    is_edge_bd = find_is_edge_bd(edges, bd_edges, size(edges,1), size(bd_edges,1));
    bd_edge_ids =find(is_edge_bd>0)';

    %% Generate result structure
    mesh = struct('nodes', vert, 'edges', double(edges), 'elements', double(tri), 'bd_edge_ids', bd_edge_ids, ...
                  'domain', domain);
end


function is_edge_bd = find_is_edge_bd(edge, bd_edge, ne, nb)
    edge_idx                = sort(edge,    2) * [ne; 1];
    bd_edge_idx             = sort(bd_edge, 2) * [ne; 1];
    [~, bd_edge_idx]        = ismember(bd_edge_idx, edge_idx);
    is_edge_bd              = zeros(ne, 1);
    is_edge_bd(bd_edge_idx) = ones(nb, 1);
end