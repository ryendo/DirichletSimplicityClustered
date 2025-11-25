function tri2edge = find_tri2edge(tri,edge)
    nt = size(tri,1);
    ne = size(edge,1);
    tri2edge = zeros(nt, 3);
    edge_idx = sort(edge, 2) * [ne; 1];
    for k = 1:nt
        edge_local = [2 3; 1 3; 1 2];
        value = sort(reshape(tri(k, edge_local), 3, 2), 2) *[ne; 1];
        [~, idx] = ismember(value, edge_idx);
        tri2edge(k,:) = idx';
    end
end