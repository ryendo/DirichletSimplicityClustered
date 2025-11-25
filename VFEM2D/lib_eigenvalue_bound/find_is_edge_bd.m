function is_edge_bd = find_is_edge_bd(edge, bd_edge, ne, nb)
    edge_idx    = sort(edge,    2) * [ne; 1];
    bd_edge_idx = sort(bd_edge, 2) * [ne; 1];
    
    [~, bd_edge_idx] = ismember(bd_edge_idx, edge_idx);
    
    is_edge_bd  = zeros(ne, 1);
    is_edge_bd(bd_edge_idx) = ones(nb, 1);
end