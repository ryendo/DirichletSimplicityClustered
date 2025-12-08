function my_intlab_config

    global INTERVAL_MODE;    
    INTERVAL_MODE=0;

    %The path of INTLAB toolbox and initialization.
    addpath("/path/to/your/INTLAB")
    
    %The path   of the library of verified eigenvalue estimation for matrix.
    addpath('../veigs')
    addpath('mesh')
    addpath('VFEM2D')
    addpath('VFEM2D/lib_eigenvalue_bound')
    addpath('VFEM2D_revised')

    
    %The path of the codes for switch between verified computing and approximate computing.
    addpath('mode_swith_interface')
    
    %The path of the library of FEMs on triangular domains
    addpath('FEM_Functions')
    
    try
        startintlab;
    catch ME
        error('An error occurred while executing startintlab: %s', ME.message);
    end

end
