function recovered_x=x_recover_dc_Pk(x,bdry_list,all_dof_size)
    
    x=intval(x);
    [~,idx_size] = size(x);
    not_on_bdry_list = setdiff(1:all_dof_size,bdry_list);
    recovered_x=intval(zeros(all_dof_size,idx_size));
    recovered_x(not_on_bdry_list,:)=x;
    
end