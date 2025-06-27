function my_intlab_config_alone

    %The path of INTLAB toolbox and initialization.
    addpath("Each_Process/Intlab_Group/Intlab_V12")
    
    %The path of the library of verified eigenvalue estimation for matrix.
    addpath('Each_Process/verified_eig_estimation')
    addpath('Each_Process/HighOrderFEM_CGYOU_2016')
    addpath('Each_Process/mat_quotients')
    
    %The path of the codes for switch between verified computing and approximate computing.
    addpath('Each_Process/mode_swith_interface')
    
    %The path of the library of FEMs on triangular domains
    addpath('Each_Process/functions')
    startintlab;
    
    global INTERVAL_MODE;
    
    INTERVAL_MODE=1;

end