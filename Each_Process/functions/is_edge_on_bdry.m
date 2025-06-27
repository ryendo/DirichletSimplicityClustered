function edge_bdry_list = is_edge_on_bdry(mesh,~)
    
    [edges_n,~] = size(mesh.edges);
    edges_in_vec = [mesh.nodes(mesh.edges(:,1),:),mesh.nodes(mesh.edges(:,2),:)];

    %Find the edge index of element which is located on the specified edge (by edge_ind) of triangle.
    bd_edge = mesh.domain( mod((1:3),3)+1,:) - mesh.domain( mod((1:3)-1,3)+1,:);
    vecs1=[edges_in_vec(:,[1,3])-  repmat(mesh.domain( mod(1,3)+1,1),edges_n,1), edges_in_vec(:,[2,4])- repmat(mesh.domain( mod(1,3)+1,2),edges_n,1)];
    vecs2=[edges_in_vec(:,[1,3])-  repmat(mesh.domain( mod(2,3)+1,1),edges_n,1), edges_in_vec(:,[2,4])- repmat(mesh.domain( mod(2,3)+1,2),edges_n,1)];
    vecs3=[edges_in_vec(:,[1,3])-  repmat(mesh.domain( mod(3,3)+1,1),edges_n,1), edges_in_vec(:,[2,4])- repmat(mesh.domain( mod(3,3)+1,2),edges_n,1)];
    
    
    test1 = abs(vecs1*[-bd_edge(1,2),-bd_edge(1,2),bd_edge(1,1),bd_edge(1,1)]')<1E-6;
    test2 = abs(vecs2*[-bd_edge(2,2),-bd_edge(2,2),bd_edge(2,1),bd_edge(2,1)]')<1E-6;
    test3 = abs(vecs3*[-bd_edge(3,2),-bd_edge(3,2),bd_edge(3,1),bd_edge(3,1)]')<1E-6;
    
    test = test1 | test2 | test3;
    
    edge_bdry_list = test';
    
end