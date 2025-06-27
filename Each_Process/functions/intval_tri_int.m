function val = intval_tri_int(f, mesh)

    % Get the number of triangles and node coordinates
    tri_n = size(mesh.elements, 1);
    nodes = mesh.nodes;
    elements = mesh.elements;

    % Get the vertices of each triangle
    p1 = nodes(elements(:,1), :);
    p2 = nodes(elements(:,2), :);
    p3 = nodes(elements(:,3), :);

    % Vectorized area calculation for each triangle
    SI = abs((p1(:,1) .* (p2(:,2) - p3(:,2)) + p2(:,1) .* (p3(:,2) - p1(:,2)) + p3(:,1) .* (p1(:,2) - p2(:,2))) / 2);

    % Calculate bounding box for each triangle
    Ix_min = min([p1(:,1), p2(:,1), p3(:,1)], [], 2);
    Ix_max = max([p1(:,1), p2(:,1), p3(:,1)], [], 2);
    Iy_min = min([p1(:,2), p2(:,2), p3(:,2)], [], 2);
    Iy_max = max([p1(:,2), p2(:,2), p3(:,2)], [], 2);

    % Initialize the function integral
    fI = intval(zeros(tri_n, 1));

    % Loop to evaluate the function f within the bounding box of each triangle
    Ix = hull(Ix_min, Ix_max);
    Iy = hull(Iy_min, Iy_max);
    fI = f(Ix, Iy);

    % Calculate the final integral value
    val = sum(fI .* SI);
end
