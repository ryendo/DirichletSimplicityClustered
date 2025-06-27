function [A,B]=create_cg_tauAB(mesh,t,func_data)
    
    Txx=func_data.Txx; Tyy=func_data.Tyy;
    
    %create A and B
    [nodes_n,~] = size(mesh.nodes);
    [elements_n,~] = size(mesh.elements);
    
    A=t^2*Txx+Tyy;
    B=intval(sparse(nodes_n,nodes_n));
    local_B=intval(ones(3,3))+intval(eye(3,3));
    for idx=1:elements_n
        element_edges = mesh.nodes( mesh.elements(idx,[3,1,2]), : ) - mesh.nodes( mesh.elements(idx,[2,3,1]), :); %3 by 2
        S =  abs(element_edges(1,:)/2*[element_edges(2,2); -element_edges(2,1)]);
        global_idx = mesh.elements(idx,1:3);
        B(global_idx, global_idx) = B(global_idx, global_idx) + local_B;
    end

    B=S/intval(12)*B;
end