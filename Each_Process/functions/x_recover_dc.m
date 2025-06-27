function u_vec=x_recover_dc(x,mesh,N)
    
    [nodes_n,~]=size(mesh.nodes);

    bd_node_list = is_node_on_bdry(mesh,N);

    dof_idx=1:nodes_n;
    dof_idx(bd_node_list)=[];

    u_vec=intval(zeros(nodes_n,1));
    u_vec(dof_idx)=intval(x);

    
    % % recovered_x=x.*is_node_on_bdry_list;
    % 
    % idx=1;
    % for i=1:nodes_n
    %     if is_node_on_bdry_list(i)
    %         recovered_x(i)=intval(0);
    %     else
    %         recovered_x(i)=x(idx);
    %         idx=idx+1;
    %     end
    % end
    
end