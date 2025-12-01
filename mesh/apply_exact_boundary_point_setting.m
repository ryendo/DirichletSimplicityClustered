function apply_exact_boundary_point_setting(mesh, tri_vertices)
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
% Desciption:
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
    p1 = tri_vertices(1, :);
    p2 = tri_vertices(2, :);
    p3 = tri_vertices(3, :);

    nodes = mesh.nodes;
    bedges = mesh.boundary_edges;

    % Tolerance: distance to vertex threshold
    vertex_tol = 1e-10;
    % Tolerance to determine which triangle edge the boundary edge is on
    edge_tol = 1e-8;

    % Loop over boundary edges
    for k = 1:size(bedges,1)

        % boundary edge node indices
        k1 = bedges(k,1);
        k2 = bedges(k,2);

        x1 = nodes(k1,:);
        x2 = nodes(k2,:);

        % Determine triangle edge by minimum distance sum
        % Distances to edges: sum of distances of both endpoints to the line
        d12 = point_to_line_distance(x1,p1,p2) + point_to_line_distance(x2,p1,p2);
        d23 = point_to_line_distance(x1,p2,p3) + point_to_line_distance(x2,p2,p3);
        d31 = point_to_line_distance(x1,p3,p1) + point_to_line_distance(x2,p3,p1);

        [min_distance, edgeID] = min([d12, d23, d31]);

        if min_distance > edge_tol
            continue
        end

        % Pick correct segment endpoints
        switch edgeID
            case 1
                a = p1;  b = p2;
            case 2
                a = p2;  b = p3;
            case 3
                a = p3;  b = p1;
        end


        % Project each node onto the correct edge
        nodes(k1,:) = project_to_segment(nodes(k1,:), a, b);
        nodes(k2,:) = project_to_segment(nodes(k2,:), a, b);

        % If close to a vertex → snap exactly
        if norm(nodes(k1,:) - a) < vertex_tol, nodes(k1,:) = a; end
        if norm(nodes(k1,:) - b) < vertex_tol, nodes(k1,:) = b; end

        if norm(nodes(k2,:) - a) < vertex_tol, nodes(k2,:) = a; end
        if norm(nodes(k2,:) - b) < vertex_tol, nodes(k2,:) = b; end

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