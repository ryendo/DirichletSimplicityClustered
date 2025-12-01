function my_intlab_config_server()

    global INTERVAL_MODE;    
    INTERVAL_MODE=1;
    global mesh_path
    mesh_path = '/tmp/';

    %The path of INTLAB toolbox and initialization.
    % addpath("/path/to/your/INTLAB")
    addpath("/home/xfliu/App/Intlab_V14")
    
    %The path   of the library of verified eigenvalue estimation for matrix.
    addpath('../veigs')
    addpath('mesh')
    addpath('VFEM2D_revised')
    addpath('VFEM2D/lib_eigenvalue_bound')
    addpath('VFEM2D')
    addpath('VFEM2D/lib_mesh')
    addpath('routines')

    
    %The path of the codes for switch between verified computing and approximate computing.
    addpath('routines/mode_swith_interface')
    
    %The path of the library of FEMs on triangular domains
    addpath('routines/FEM_Functions')

    addpath('EigenvaluesInterpolation')
    
    try
        startintlab;
    catch ME
        error('An error occurred while executing startintlab: %s', ME.message);
    end

end
