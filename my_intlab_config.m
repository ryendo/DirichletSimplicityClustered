function my_intlab_config

    global INTERVAL_MODE;    
    INTERVAL_MODE=1;

    %The path of INTLAB toolbox and initialization.
    % addpath("/path/to/your/INTLAB")
    
    %The path of the library of verified eigenvalue estimation for matrix.
    addpath('veigs')
    addpath('VFEM2D_revised')
    addpath('routines')

    
    %The path of the codes for switch between verified computing and approximate computing.
    addpath('mode_swith_interface')
    
    %The path of the library of FEMs on triangular domains
<<<<<<< Updated upstream:my_intlab_config_alone.m
    addpath('routines/FEM_Functions')
=======
    addpath('FEM_Functions')
>>>>>>> Stashed changes:my_intlab_config.m
    
    try
        startintlab;
    catch ME
        error('An error occurred while executing startintlab: %s', ME.message);
    end

end
