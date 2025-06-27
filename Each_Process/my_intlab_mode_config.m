function my_intlab_mode_config(n)

    %The path of INTLAB toolbox and initialization.
    addpath("Intlab_Group/Intlab_V12_no"+num2str(n))
    
    %The path of the library of verified eigenvalue estimation for matrix.
    addpath('verified_eig_estimation')
    addpath('HighOrderFEM_CGYOU_2016')
    addpath('mat_quotients')
    
    %The path of the codes for switch between verified computing and approximate computing.
    addpath('mode_swith_interface')
    
    %The path of the library of FEMs on triangular domains
    addpath('functions')
    startintlab;
    
    global INTERVAL_MODE;
    
    INTERVAL_MODE=1;

end