function blams = calc_eigen_bounds(tri_intval,N)
    x=tri_intval(5); y=tri_intval(6);
    
    if x<0.7
        
        mesh=get_mesh_for_cg(tri_intval,N);
        
        
        func_data=init_func(mesh,N);
        
    else
        
        mesh=get_mesh_for_cg_right_triangles(tri_intval,N);
        
        
        func_data=init_func(mesh,N);
        
    end
    
    [A_cg_eigs,B_cg_eigs]=create_cg_AB(mesh,func_data);
    
    [A_cg_dc_no_bdry,B_cg_dc_no_bdry]=add_dc_cg(A_cg_eigs,B_cg_eigs,false,mesh,N);
    
    cg_lambdas=veig(A_cg_dc_no_bdry,B_cg_dc_no_bdry,1:4);
    
    % Compute upper and lower bounds
    Ch = intval(0.493)/N;
    ulams=cg_lambdas;
    llams=(cg_lambdas./(1+Ch^2.*cg_lambdas));
    lams=hull(ulams,llams);
    blams=lams;

end