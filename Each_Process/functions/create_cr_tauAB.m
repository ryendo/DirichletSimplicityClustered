function [A,B]=create_cr_tauAB(mesh,t,~)
    [edges_n,~] = size(mesh.edges);
    [elements_n,~] = size(mesh.elements);
    
    global Txxcr; global Tyycr;
    %create the table of element-mesh.edges relation.
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
    
    % create A and B
    A=t^2*Txxcr+Tyycr;
    B=intval(sparse(edges_n,edges_n));
    for idx=1:elements_n
        element_edges = mesh.nodes( mesh.elements(idx,[3,1,2]), : ) - mesh.nodes( mesh.elements(idx,[2,3,1]), :); %3 by 2
        S =  abs(element_edges(1,:)/2*[element_edges(2,2); -element_edges(2,1)]);
        
        % contribute local matirix to A and B
        edge_index = tri_by_edge(idx,:);
        B(edge_index, edge_index) = B(edge_index, edge_index) + S*intval(eye(3,3))/3;
    end
end

