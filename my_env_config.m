function my_env_config()

    global INTERVAL_MODE;    
    INTERVAL_MODE=1;

    global gmsh_command 
    gmsh_command = '/opt/homebrew/bin/gmsh'  %Path of gmsh command.;
    %The computation was tested with gmsh version 4.8.4.

    global mesh_path
    mesh_path = '/tmp/'     %Place to save temporary mesh files.


    %The path of INTLAB toolbox and initialization.
    %addpath("/path/to/your/INTLAB")
    addpath("/User/xfliu/App/Intlab_V14")
    
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
