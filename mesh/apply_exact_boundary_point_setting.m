function mesh=apply_exact_boundary_point_setting(mesh, tri_vertices)
% apply_exact_boundary_point_setting - Adjust boundary-node coordinates so
% that all boundary nodes lie exactly on the triangle boundary.
%
% INPUT:
%   mesh           : a structure containing fields:
%                       mesh.nodes            (N × 2)
%                       mesh.boundary_edges   (M × 2)
%   tri_vertices   : 3 × 2 matrix, each row is (x,y) of triangle vertex
%
% NOTE:
%   The function modifies mesh.nodes in-place if mesh is a handle-like struct,
%   otherwise return the modified mesh.
%
% Description:
% Given 'mesh' is a triangulation of domain specified by 3*2 tri_vertices,
% each row of which is a vertex.
% This function proccess the boundary nodes in mesh.nodes to make sure the
% boundary nodes is really on the border of domain. 
% Note that mesh.boundary_edges provides the list of boundary edges, each
% edge is represented by the indices of its two end vertices.

% For example, let p1, p2, p3 be the 3 vertices specified by tri_vertices.
% For an edge [k1,k2], use the direction to dertermined the domain edge on which 
% edge [k1, k2] is locaed on. Suppose nodes(k1,:) and nodes(k2,:) is on triangle domain
% edge p1-p2, then for each (x,y) = nodes(k1,:) and (x,y)=nodes (k2,:), y coordinate
% is recalculated by the line determined by p1-p2.
% If nodes(k1,:) is very near to vertex p of p1, p2 or p3, then set the value of
% nodes(k1,:) by p.


    % Triangle vertices
    p1_ori = tri_vertices(1, :);
    p2_ori = tri_vertices(2, :);
    p3_ori = tri_vertices(3, :);
    p1 = I_mid(p1_ori);
    p2 = I_mid(p2_ori);
    p3 = I_mid(p3_ori);

    nodes = I_intval(mesh.nodes);
    bedges = mesh.boundary_edges;

    % Tolerance: distance to vertex threshold
    vertex_tol = 1e-6;
    % Tolerance to determine which triangle edge the boundary edge is on
    edge_tol = 1e-6;

    % Loop over boundary edges
    x_idx_list = unique(bedges(:))';
    for idx = x_idx_list
        x = I_mid(nodes(idx,:));

        % If close to a vertex → snap exactly
        if norm(nodes(idx,:) - p1) < vertex_tol
	    nodes(idx,:) = p1_ori; continue;
	end
        if norm(nodes(idx,:) - p2) < vertex_tol
	    nodes(idx,:) = p2_ori; continue;
	end
        if norm(nodes(idx,:) - p3) < vertex_tol
	    nodes(idx,:) = p3_ori; continue;
	end

        % Determine triangle edge by minimum distance
        % Distances to edges
        d12 = point_to_line_distance(x,p1,p2);
        d23 = point_to_line_distance(x,p2,p3);
        d31 = point_to_line_distance(x,p3,p1);

        [min_distance, edgeID] = min([d12, d23, d31]);

        if min_distance > edge_tol
            continue
        end

        % Pick correct segment endpoints
        switch edgeID
            case 1
                a = p1_ori;  b = p2_ori;
            case 2
                a = p2_ori;  b = p3_ori;
            case 3
                a = p3_ori;  b = p1_ori;
        end

        % Project each node onto the correct edge
        nodes(idx,:) = project_to_segment(nodes(idx,:), a, b);

    end

    % Write back
    mesh.nodes = nodes;

end


%% ===============================================================
%   Helper Functions
% ===============================================================

% Distance from point P to infinite line AB
function d = point_to_line_distance(P, A, B)
    AP = P - A;
    AB = B - A;
    d = norm(AP - (dot(AP,AB)/dot(AB,AB))*AB);
end

% Project point P onto segment AB (clamped projection)
function Pproj = project_to_segment(P, A, B)
    AB = B - A;
    t = dot(P-A, AB) / dot(AB, AB);
    % clamp to [0,1]
    t = max(0, min(1, t));
    Pproj = A + t * AB;
end