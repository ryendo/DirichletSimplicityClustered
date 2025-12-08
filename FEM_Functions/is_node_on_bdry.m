function edge_bdry_list = is_node_on_bdry(mesh,~)
    
    [nodes_n,~] = size(mesh.nodes);
    %Find the edge index of element which is located on the specified edge (by edge_ind) of triangle.
    edge_ind = 1:3;
    node_ind = 1:nodes_n;
    edge = mesh.domain( mod(edge_ind+1,3)+1,:) - mesh.domain( mod(edge_ind,3)+1,:);

    [node_ind_n,~] = size(node_ind);
    vecs1=[mesh.nodes(node_ind,1) - repmat(mesh.domain( mod(1,3)+1,1),node_ind_n,1), mesh.nodes(node_ind,2) - repmat(mesh.domain( mod(1,3)+1,2),node_ind_n,1)];
    vecs2=[mesh.nodes(node_ind,1) - repmat(mesh.domain( mod(2,3)+1,1),node_ind_n,1), mesh.nodes(node_ind,2) - repmat(mesh.domain( mod(2,3)+1,2),node_ind_n,1)];
    vecs3=[mesh.nodes(node_ind,1) - repmat(mesh.domain( mod(3,3)+1,1),node_ind_n,1), mesh.nodes(node_ind,2) - repmat(mesh.domain( mod(3,3)+1,2),node_ind_n,1)];
    
    test1 = abs(vecs1*[-edge(1,2),edge(1,1)]')<1E-6;
    test2 = abs(vecs2*[-edge(2,2),edge(2,1)]')<1E-6;
    test3 = abs(vecs3*[-edge(3,2),edge(3,1)]')<1E-6;
    test = test1 | test2 | test3;
    
    edge_bdry_list = test';
    
end

