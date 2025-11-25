function my_intlab_config_alone

    %The path of INTLAB toolbox and initialization.
    % addpath("/path/to/your/INTLAB")
    addpath("/home/rendo/Documents/Intlab_V12")
    
    %The path of the library of verified eigenvalue estimation for matrix.
    addpath('veigs')
    addpath('VFEM2D_revised')
    addpath('Each_Process')
    addpath('Each_Process/mat_quotients')
    
    %The path of the codes for switch between verified computing and approximate computing.
    addpath('Each_Process/mode_swith_interface')
    
    %The path of the library of FEMs on triangular domains
    addpath('Each_Process/FEM_Functions')
    startintlab;
    
    global INTERVAL_MODE;
    
    INTERVAL_MODE=1;

end
