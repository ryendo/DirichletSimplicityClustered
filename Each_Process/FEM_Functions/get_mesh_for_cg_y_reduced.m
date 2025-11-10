function mesh = get_mesh_for_cg_y_reduced(tri_node, N_h, N_v)
    % Generates a structured mesh for a given triangle.
    %
    % This function creates a mesh by dividing the triangle into horizontal strips
    % and populating nodes on each level, then connecting them to form elements.
    % This corrected version avoids node duplication and ensures a connected mesh.
    % It is also compatible with intval inputs.
    %
    % INPUTS:
    %   tri_node: 3x2 matrix (double or intval) with coordinates of the triangle vertices.
    %   N_h:      Number of divisions along the base of the triangle.
    %   N_v:      Number of vertical divisions.
    
    if nargin < 3
        N_v = N_h;
    end

    % Find the top vertex index, compatible with intval objects.
    % For intval, we cannot get the index from max() directly.
    % Instead, we find the index of the maximum of the midpoints of the intervals.
    y_coords_mid = I_mid(tri_node(:,2));
    [~, top_idx] = max(y_coords_mid);
    
    other_indices = setdiff(1:3, top_idx);
    p3 = tri_node(top_idx, :);
    p1 = tri_node(other_indices(1), :);
    p2 = tri_node(other_indices(2), :);

    % Ensure p1 is to the left of p2 for consistent ordering, using midpoints
    if I_mid(p1(1)) > I_mid(p2(1))
        [p1, p2] = deal(p2, p1); % Swap p1 and p2
    end
    
    x1 = p1(1); y1 = p1(2);
    x2 = p2(1); y2 = p2(2);
    x3 = p3(1); y3 = p3(2);

    % Pre-allocate memory
    total_nodes = (N_h + 2) * (N_v + 1); % Safe over-estimation
    nodes = I_intval(zeros(round(total_nodes), 2)); % Use your intval constructor
    row_indices = cell(N_v + 1, 1);
    node_counter = 0;

    % --- Unified Node Generation ---
    for r = 0:N_v
        % y-coordinate for the current row (from top y3 to bottom average y)
        y_level = y3 - (y3 - (y1+y2)/2) * (r / N_v);
        
        % Determine x_start and x_end for this y_level by linear interpolation
        % along the left (p3-p1) and right (p3-p2) edges.
        if I_mid(y3) == I_mid(y1)
             x_start = x1;
        else
            x_start = x3 + (x1 - x3) * (y3 - y_level) / (y3 - y1);
        end
        
        if I_mid(y3) == I_mid(y2)
            x_end = x2;
        else
            x_end = x3 + (x2 - x3) * (y3 - y_level) / (y3 - y2);
        end

        % Number of nodes in this row
        num_row_nodes = 1 + round(r * N_h / N_v);
        
        x_vals = linspace(x_start, x_end, num_row_nodes);
        
        current_row_node_ids = zeros(1, num_row_nodes);
        for i = 1:num_row_nodes
            node_counter = node_counter + 1;
            nodes(node_counter, :) = [x_vals(i), y_level];
            current_row_node_ids(i) = node_counter;
        end
        row_indices{r + 1} = current_row_node_ids;
    end
    nodes = nodes(1:node_counter, :); % Trim unused pre-allocated space

    % --- Unified Triangulation (Element Generation) ---
    elements = [];
    for r = 1:N_v
        P = row_indices{r};   % Top row of the strip
        Q = row_indices{r+1}; % Bottom row of the strip
        
        i = 1; % Pointer for P
        j = 1; % Pointer for Q
        
        while i < length(P) || j < length(Q)
            is_p_end = (i == length(P));
            is_q_end = (j == length(Q));
            
            if is_p_end
                elements = [elements; P(i), Q(j), Q(j+1)];
                j = j + 1;
            elseif is_q_end
                elements = [elements; P(i), Q(j), P(i+1)];
                i = i + 1;
            else
                if I_mid(nodes(P(i+1), 1)) <= I_mid(nodes(Q(j+1), 1))
                    elements = [elements; P(i), Q(j), P(i+1)];
                    i = i + 1;
                else
                    elements = [elements; P(i), Q(j), Q(j+1)];
                    j = j + 1;
                end
            end
        end
    end

    % --- Create edges from element connectivity ---
    edge_list = [elements(:,[1,2]); elements(:,[2,3]); elements(:,[3,1])];
    edge_list = sort(edge_list, 2);
    edges = unique(edge_list, 'rows');
    
    % --- Define domain boundary ---
    domain = [p1; p2; p3];

    mesh = struct('nodes', nodes, 'edges', edges, 'elements', elements, 'domain', domain);
end