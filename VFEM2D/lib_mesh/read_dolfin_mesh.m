function mesh = read_dolfin_mesh(filename)
% READ_DOLFIN_MESH Generate a mesh structure from a loaded Dolfin XML mesh file
%
%   mesh = read_dolfin_mesh(filename)
%
% Inputs:
%   filename - Path to the XML file to be read (e.g., 'mesh.xml')
%
% Outputs:
%   mesh - Structure (fields: nodes, edges, elements, domain, min_edge_length)
%          nodes: Each row contains [x, y] coordinates (z component in XML is ignored)
%          elements: Each row contains 3 vertex indices (XML triangle elements; Note: XML is 0-based, so +1 is added for MATLAB)
%          edges: Unique edges extracted from elements (each edge sorted in ascending order)
%          domain: Stores vertices of the first triangle in [1, 3, 2] order as a sample
%          min_edge_length: The minimum value among the maximum edge lengths of each triangle element
%
% Example:
%   mesh = read_dolfin_mesh('mesh.xml');
%  Revision history:
%   2024/06/10: Xuefeng Liu, First edition by ChatGPT based on user request

    % Read XML file
    xDoc = xmlread(filename);
    
    %% Extract vertex information
    verticesList = xDoc.getElementsByTagName('vertex');
    numVertices = verticesList.getLength;
    nodes = zeros(numVertices, 2);
    
    for k = 0:numVertices-1
        vertex = verticesList.item(k);
        % Get x, y attribute values as strings and convert to numbers
        xStr = char(vertex.getAttribute('x'));
        yStr = char(vertex.getAttribute('y'));
        nodes(k+1, :) = [str2double(xStr), str2double(yStr)];
    end
    
    %% Extract triangle (cell) information
    trianglesList = xDoc.getElementsByTagName('triangle');
    numTriangles = trianglesList.getLength;
    elements = zeros(numTriangles, 3);
    
    for k = 0:numTriangles-1
        triangle = trianglesList.item(k);
        v0 = str2double(char(triangle.getAttribute('v0')));
        v1 = str2double(char(triangle.getAttribute('v1')));
        v2 = str2double(char(triangle.getAttribute('v2')));
        % Indices in XML are 0-based, so +1 for MATLAB
        elements(k+1, :) = [v0+1, v1+1, v2+1];
    end
    
    %% Create edge information (extract edges from each triangle and make them unique)
    edge_list = [];
    for k = 1:size(elements,1)
        tri = elements(k,:);
        % Each triangle has 3 edges: (v1,v2), (v2,v3), (v3,v1)
        edges_tri = [tri([1,2]); tri([2,3]); tri([3,1])];
        % Sort each edge so the smaller index comes first
        edges_tri = sort(edges_tri,2);
        edge_list = [edge_list; edges_tri];
    end
    edges = unique(edge_list, 'rows');
    
    %% Calculate the maximum edge length for each triangle element, then find the minimum among them
    triangle_max_lengths = zeros(numTriangles, 1);
    for k = 1:numTriangles
        tri = elements(k,:);
        % Calculate length of each edge
        d1 = norm(nodes(tri(1),:) - nodes(tri(2),:));
        d2 = norm(nodes(tri(2),:) - nodes(tri(3),:));
        d3 = norm(nodes(tri(3),:) - nodes(tri(1),:));
        triangle_max_lengths(k) = max([d1, d2, d3]);
    end
    min_edge_length = min(triangle_max_lengths);
    
    %% Define domain
    % As an example here, the domain is set using the first triangle's vertices in [1, 3, 2] order
    %    domain = nodes([1, 3, 2], :);
    domain = [];

    bd_edge_ids = get_boundary_edge(elements, edges);
    
    %% Generate result structure
    mesh = struct('nodes', nodes, 'edges', double(edges), 'elements', double(elements), ...
                  'bd_edge_ids', bd_edge_ids, ... 
                  'domain', domain, 'min_edge_length', min_edge_length);
end

function boundary_edge_ids = get_boundary_edge(elements, edges)
    % GET_BOUNDARY_EDGE_INDICES
    %   boundary_edge_ids = get_boundary_edge_indices(mesh)
    %
    % Inputs:
    %   mesh - struct (fields: nodes, edges, elements, domain)
    %          elements: [M x 3] matrix (vertex indices of triangles, 1-based)
    %          edges: [E x 2] matrix (unique edges, each edge sorted in ascending order)
    %
    % Outputs:
    %   boundary_edge_ids - List of indices in mesh.edges corresponding to boundary edges
    %
    % Algorithm:
    %   1. Extract 3 edges from each triangle to create the full edge list 'allEdges'.
    %   2. Determine the occurrence count of each unique edge using unique and accumarray.
    %   3. Identify edges with an occurrence count of 1 as boundary edges.
    %   4. Check which edges in mesh.edges correspond to the boundary edges and return their indices.

    % --- 1. Enumerate edges for each triangle ---
    numElements = size(elements, 1);
    allEdges = zeros(numElements * 3, 2);
    idx = 1;
    for i = 1:numElements
        tri = elements(i,:);
        % Edges of each triangle (sorted in ascending order)
        edge1 = sort([tri(1), tri(2)]);
        edge2 = sort([tri(2), tri(3)]);
        edge3 = sort([tri(3), tri(1)]);
        allEdges(idx, :)   = edge1;
        allEdges(idx+1, :) = edge2;
        allEdges(idx+2, :) = edge3;
        idx = idx + 3;
    end

    % --- 2. Get unique edges and their occurrence counts ---
    [uEdges, ~, ic] = unique(allEdges, 'rows');
    counts = accumarray(ic, 1);
    
    % --- 3. Edges with an occurrence count of 1 are boundary edges ---
    boundaryEdges = uEdges(counts == 1, :);
    
    % --- 4. Extract indices within mesh.edges corresponding to boundary edges ---
    numMeshEdges = size(edges, 1);
    boundary_edge_ids = [];
    for j = 1:numMeshEdges
        % Assuming mesh.edges(j,:) is already sorted in ascending order
        current_edge = edges(j,:);
        if ismember(current_edge, boundaryEdges, 'rows')
            boundary_edge_ids(end+1,1) = j; %#ok<AGROW>
        end
    end
end