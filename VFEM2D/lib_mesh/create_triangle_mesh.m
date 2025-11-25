function mesh = create_triangle_mesh(varargin)
% CREATE_TRIANGLE_MESH Generates a uniform triangular mesh by subdividing a given triangle.
%
%   mesh = create_triangle_mesh(tri_node, n)
%
% Description:
%   This function subdivides a single large triangle into
%   a uniform mesh of smaller triangles. The subdivision strategy 
%   depends on the geometry of the input triangle (specifically the largest angle) 
%   to ensure the quality of the resulting elements, but the grid spacing 
%   remains uniform based on the subdivision parameter 'n'.
%
%   Note: This code utilizes Interval Arithmetic functions (I_intval, I_sup, I_pi).
%
% Inputs:
%   varargin{1} (tri_node): [3 x 2] matrix representing the (x,y) coordinates
%                           of the vertices of the initial triangle.
%   varargin{2} (n)       : Integer. The subdivision parameter. The triangle
%                           edges will be divided into 'n' uniform segments.
%
% Outputs:
%   mesh: A structure containing the generated mesh data with fields:
%         - nodes: [N x 2] coordinates of mesh vertices.
%         - edges: [E x 2] indices of graph edges.
%         - elements: [M x 3] vertex indices for each triangle element.
%         - bd_edge_ids: Indices of edges located on the boundary.
%         - domain: The coordinates of the original triangle (reordered).
%         - min_edge_length: The minimum of the maximum edge lengths of the elements.

    tri_node = varargin{1};
     
    % --- 1. Geometric Analysis & Reordering ---
    % Calculate edge vectors and lengths of the input triangle
    edges = tri_node([3,1,2],:) - tri_node([2,3,1],:);
    lens = sqrt( edges(:,1).^2 + edges(:,2).^2);
    
    % Find the longest edge and its corresponding index
    % I_sup is used here to get the upper bound of the interval for comparison
    [~, ind] = max(I_sup(lens));
    ind = ind(1,1); % 'ind' is the index of the vertex opposite to the maximal edge.
    
    % Reorder nodes so that the reference vertex (ind) corresponds to the
    % angle opposite the longest edge.
    node_ind = [ind, mod(ind,3)+1, mod(ind+1,3)+1];
    domain = tri_node(node_ind, :);
    
    % Recalculate edges based on the reordered domain
    edges = domain([3,1,2],:) - domain([2,3,1],:);

    % Calculate the cosine of the angle at the reference node using Law of Cosines
    % cos(A) = (b^2 + c^2 - a^2) / 2bc
    tmp = (sum(lens(node_ind(2:3)).^2) - lens(ind)^2) / I_intval(2) / lens(node_ind(2)) / lens(node_ind(3));
    max_angle = acos(tmp);
       
    n = varargin{2};
    
    % Define step sizes (h1, h2) for grid generation using interval arithmetic
    h1 = I_intval(1)/n; % Step size along e3 direction (side adjacent to ref node)
    h2 = I_intval(1)/n; % Step size along e2 direction (other side adjacent to ref node)
    vec1 = edges(3,:);  % Vector along one side
    vec2 = -edges(2,:); % Vector along the other side

    elements = [];
    edges_list = []; % Renamed from 'edges' to avoid conflict with previous variable
    start = 0;
    ref_p = tri_node(ind, :);

    % --- 2. Mesh Generation ---
    % The strategy changes based on whether the triangle is acute/right or obtuse.

    % Case A: Angle is less than 90 degrees (Acute or Right triangle strategy)
    % Uses standard uniform subdivision.
    if max_angle < I_pi/2
    
        nodes = [];
        % Generate Grid Nodes
        for k = 0:n
            x = h1 * (0:n-k);
            shift_p = (ref_p + k*h2*vec2);
            tmp = x' * vec1;
            new_nodes = [tmp(:,1)+shift_p(1), tmp(:,2)+shift_p(2)];
            nodes = [nodes; new_nodes];
        end
    
        % Generate Connectivity (Elements and Edges)
        for k = 1:n
            % Indices offsets for the current row
            basic_tri1 = [1, 2, n+2-(k-1)];
            basic_tri2 = [2, n+2-(k-1)+1, n+2-(k-1)];
            
            % Create "upward" pointing triangles
            for l = 1:n-(k-1)
               edges_list = [edges_list; basic_tri1(1:2)+l-1+start]; 
               edges_list = [edges_list; basic_tri1(2:3)+l-1+start]; 
               edges_list = [edges_list; basic_tri1([1,3])+l-1+start]; 

               elements = [elements; basic_tri1+l-1+start];
            end
            
            % Create "downward" pointing triangles (filling the gaps)
            for l = 1:n-1-(k-1)
                elements = [elements; basic_tri2+l-1+start];
            end
            
            % Update start index for the next row of nodes
            start = start + (n+1) - (k-1);
        end

    % Case B: Angle is >= 90 degrees (Obtuse triangle strategy)
    % Uses a refined subdivision that adds extra nodes on the longest edge/hypotenuse
    % to improve element quality (preventing bad aspect ratios).
    else
    
        nodes = [];
        % Generate Nodes (Standard grid + Extra points)
        for k = 0:n
            x = h1 * (0:n-k);
            shift_p = (ref_p + k*h2*vec2);
            tmp = x' * vec1;
            new_nodes = [tmp(:,1)+shift_p(1), tmp(:,2)+shift_p(2)];
            nodes = [nodes; new_nodes];
            
            % Add intermediate points on the diagonal (max edge direction)
            if k < n 
                new_point_on_max_edge = (n-k-0.5)*h1*vec1 + h2*vec2*(k+0.5) + ref_p;
                nodes = [nodes; new_point_on_max_edge];
            end
        end

        % Generate Connectivity for the obtuse refinement pattern
        for k = 1:n
            
            % Define topology templates for this specific refinement
            basic_tri1 = [1, n+2-(k-1)+2, n+2-(k-1)+1];
            basic_tri2 = [1, 2, n+2-(k-1)+2];
            
            for l = 1:n-1-(k-1)
                elements = [elements; basic_tri1+l-1+start];
                elements = [elements; basic_tri2+l-1+start];
                
                edges_list = [edges_list; [1,2] + l-1+start]; 
                edges_list = [edges_list; [1,n+2-(k-1)+2] + l-1+start]; 
                edges_list = [edges_list; [1,n+2-(k-1)+1] + l-1+start]; 
            end
            
            % Handle the closure of the row/diagonal specific to this pattern
            tmp1 = [n-(k-1), n+2-(k-1), 2*n-2*k+4] + start;
            tmp2 = [n-(k-1), n-(k-1)+1, n+2-(k-1)] + start;            
            
            edges_list = [edges_list; tmp2([1,2])]; 
            edges_list = [edges_list; tmp2([2,3])]; 
            edges_list = [edges_list; tmp2([1,3])]; 
            edges_list = [edges_list; tmp1([2,3])]; 
            edges_list = [edges_list; tmp1([1,3])]; 
            
            elements = [elements; tmp1; tmp2];
            start = start + (n+1) - (k-1) + 1;
        end
    end
    
    % --- 3. Post-Processing: Metrics and Boundary Detection ---
    
    % Calculate the size of elements (missing in original code, calculated here)
    numTriangles = size(elements, 1); 
    triangle_max_lengths = zeros(numTriangles, 1);
    
    % Compute the maximum edge length for every generated triangle
    for k = 1:numTriangles
        tri = elements(k,:);
        d1 = norm(nodes(tri(1),:) - nodes(tri(2),:));
        d2 = norm(nodes(tri(2),:) - nodes(tri(3),:));
        d3 = norm(nodes(tri(3),:) - nodes(tri(1),:));
        triangle_max_lengths(k) = max([d1, d2, d3]);
    end
    % Identify the minimum of these maximum lengths (mesh quality metric)
    min_edge_length = min(triangle_max_lengths);

    % Identify boundary edges
    bd_edge_ids = get_boundary_edge(elements, edges_list);

    % Construct the final mesh structure
    mesh = struct('nodes', nodes, ...
                  'edges', double(edges_list), ...
                  'elements', double(elements), ...
                  'bd_edge_ids', bd_edge_ids, ... 
                  'domain', domain, ...
                  'min_edge_length', min_edge_length);
end


function boundary_edge_ids = get_boundary_edge(elements, edges)
% GET_BOUNDARY_EDGE_INDICES Extracts indices of edges that lie on the boundary.
%
% Inputs:
%   elements: [M x 3] Matrix of triangle vertex indices.
%   edges:    [E x 2] Matrix of all edges (global indexing).
%
% Outputs:
%   boundary_edge_ids: Indices in the 'edges' input array that correspond
%                      to the boundary of the mesh.

    % --- 1. Enumerate edges for each triangle ---
    numElements = size(elements, 1);
    allEdges = zeros(numElements * 3, 2);
    idx = 1;
    for i = 1:numElements
        tri = elements(i,:);
        % Edges of each triangle (sorted in ascending order to handle orientation)
        edge1 = sort([tri(1), tri(2)]);
        edge2 = sort([tri(2), tri(3)]);
        edge3 = sort([tri(3), tri(1)]);
        allEdges(idx, :)   = edge1;
        allEdges(idx+1, :) = edge2;
        allEdges(idx+2, :) = edge3;
        idx = idx + 3;
    end

    % --- 2. Get unique edges and their occurrence counts ---
    % A boundary edge appears exactly once in 'allEdges'. 
    % Internal edges appear exactly twice (shared by two triangles).
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