function func_data = init_func2(mesh,N)
    
    func_data = init_u_u_(mesh,N);

end


function func_data = init_u_u_(mesh,N)
    tic
    [nodes_n,~] = size(mesh.nodes);
    [edges_n,~] = size(mesh.edges);
    [elements_n,~] = size(mesh.elements);

    tri_by_edge = zeros(elements_n,3);
    edge_value_index = sort(mesh.edges,2)*[edges_n;1];
    for k=1:elements_n  %create the table of element-mesh.edges relation.
        for l = 1:3  %the l-th edge.
            node_ind =[ mod(l,3)+1, mod(l+1,3)+1];
            value = sort(mesh.elements(k, node_ind) )*[edges_n;1];
            [r,ind] = ismember( value, edge_value_index);
            tri_by_edge(k,l) = ind; %The direction of an edge is not counted, which is needed for other mesh.elements, e.g., Fujino-Morley FEM
        end
    end
    toc

    tic
    Txx_u=sparse(nodes_n,nodes_n); Txx_l=sparse(nodes_n,nodes_n);
    Txy_u=sparse(nodes_n,nodes_n); Txy_l=sparse(nodes_n,nodes_n);
    Tyy_u=sparse(nodes_n,nodes_n); Tyy_l=sparse(nodes_n,nodes_n);
    T_u=sparse(nodes_n,nodes_n);   T_l=sparse(nodes_n,nodes_n);
    
    for idx=1:elements_n
        element_edges = mesh.nodes( mesh.elements(idx,[3,1,2]), : ) - mesh.nodes( mesh.elements(idx,[2,3,1]), :); %3 by 2
        S =  abs(element_edges(1,:)/2*[element_edges(2,2); -element_edges(2,1)]);
        global_idx = mesh.elements(idx,1:3);
        edge_index = tri_by_edge(idx,:);
        
        e_e_xx = element_edges(:,2)*element_edges(:,2)';
        e_e_xy = element_edges(:,1)*element_edges(:,2)';
        e_e_yy = element_edges(:,1)*element_edges(:,1)';
        
        feature('setround', Inf);
        Txx_u(global_idx, global_idx)=Txx_u(global_idx, global_idx)+sup(e_e_xx/(4*S));
        Txy_u(global_idx, global_idx)=Txy_u(global_idx, global_idx)+sup(e_e_xy/(4*S));
        Tyy_u(global_idx, global_idx)=Tyy_u(global_idx, global_idx)+sup(e_e_yy/(4*S));
        T_u(global_idx, global_idx)  =T_u(global_idx, global_idx)+sup(intval((ones(3,3)+eye(3,3)))*S/intval(12));

        feature('setround', -Inf);
        Txx_l(global_idx, global_idx)=Txx_l(global_idx, global_idx)+inf(e_e_xx/(4*S));
        Txy_l(global_idx, global_idx)=Txy_l(global_idx, global_idx)+inf(e_e_xy/(4*S));
        Tyy_l(global_idx, global_idx)=Tyy_l(global_idx, global_idx)+inf(e_e_yy/(4*S));
        T_l(global_idx, global_idx)  =T_l(global_idx, global_idx)+inf(intval((ones(3,3)+eye(3,3)))*S/intval(12));

        feature('setround', 0.5);

        % display progress
        % fprintf(' basic matrices %4.2f percent prepared\n',round(idx/elements_n*100,4))
    end
    toc
    % global Txx; global Txy; global Tyy; global T; global Txxcr; global Txycr; global Tyycr; global Tcr;
    Txx=infsup(Txx_l,Txx_u); Txy=infsup(Txy_l,Txy_u); Tyy=infsup(Tyy_l,Tyy_u); T=infsup(T_l,T_u);

    func_data = struct('Txx',Txx,'Txy',Txy,'Tyy',Tyy,'T',T);
end