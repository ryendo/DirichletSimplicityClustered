function my_intlab_mode_config(n)

    % --- Automatically determine absolute paths ---
    % Get the full path to this script file ('my_intlab_mode_config.m')
    this_script_path = mfilename('fullpath');
    
    % Get the directory containing this script ('.../Each_Process')
    each_process_dir = fileparts(this_script_path);
    
    % Get the main project directory ('.../DirichletSimplicityClustered')
    project_root_dir = fileparts(each_process_dir);
    
    
    % --- Add the specific INTLAB path for this process ---
    % Dynamically build the path to the correct Intlab_V12_no<n> folder
    intlab_folder_for_this_process = "Intlab_V12_no" + n;
    addpath('/home/rendo/Documents/DirichletSimplicityClustered/Each_Process/Intlab_Group/Intlab_V12');
    
    
    % --- Add paths to other libraries using absolute paths ---
    addpath(fullfile(project_root_dir, 'verified_eig_estimation'));
    addpath(fullfile(project_root_dir, 'HighOrderFEM_CGYOU_2016_Dirichlet'));
    addpath(fullfile(each_process_dir, 'mat_quotients'));
    addpath(fullfile(each_process_dir, 'mode_swith_interface'));
    addpath(fullfile(each_process_dir, 'FEM_Functions'));
    
    % Now that the path is correct, startintlab should work
    startintlab;
    
    global INTERVAL_MODE;
    INTERVAL_MODE=1;

end
