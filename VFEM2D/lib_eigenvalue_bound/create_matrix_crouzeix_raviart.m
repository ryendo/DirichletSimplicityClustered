function [A0, A1] = create_matrix_crouzeix_raviart(tri, edge, node, tri_by_edge)
%% Function to calcuate the inner product of P1 bases of Crouzeix-Raviart FEM
%% Xuefeng LIU (2014/06/16)

    t_num = size(tri,1);
    e_num = size(edge,1); num = e_num;

    A0 = I_sparse( num, num );  
    A1 = I_sparse( num, num );


    for k=1:t_num
        t=tri(k,:);
        edges = node( t([3,1,2]), : ) - node( t([2,3,1]), : ); %3 by 2
        edge_index = tri_by_edge(k,:);        
        S = abs(0.5* edges(1,:)*[edges(2,2); -edges(2,1)]);
        e_e = edges*edges';
        A0( edge_index, edge_index ) = A0( edge_index, edge_index ) + S*eye(3,3)/3.0;
        A1( edge_index, edge_index ) = A1( edge_index, edge_index ) + e_e/S;
    end
end
