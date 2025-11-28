function my_intlab_config_alone

    global INTERVAL_MODE;    
    INTERVAL_MODE=0;

    %The path of INTLAB toolbox and initialization.
    % addpath("/path/to/your/INTLAB")
    addpath("/Users/ryoki/Documents/GitHub/ConjectureIsWrong/Intlab_V12")
    
    %The path of the library of verified eigenvalue estimation for matrix.
    addpath('veigs')
    addpath('VFEM2D_revised')
    addpath('routines')
    
    %The path of the codes for switch between verified computing and approximate computing.
    addpath('routines/mode_swith_interface')
    
    %The path of the library of FEMs on triangular domains
    addpath('routines/FEM_Functions')
    
    try
        startintlab;
    catch ME
        error('An error occurred while executing startintlab: %s', ME.message);
    end

end
