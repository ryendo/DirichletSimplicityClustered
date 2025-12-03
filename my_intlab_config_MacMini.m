global INTERVAL_MODE   
INTERVAL_MODE=1
global mesh_path
mesh_path = '/tmp/'

global region_bound_validation_file_OK region_bound_validation_file_NG

if INTERVAL_MODE > 0
    region_bound_validation_file_OK ='./results/region_0.026_to_0.1_OK_verified.txt'
    region_bound_validation_file_NG ='./results/region_0.026_to_0.1_NG_verified.txt'
else
    region_bound_validation_file_OK ='./results/region_OK.txt'
    region_bound_validation_file_NG ='./results/region_NG.txt'
end

%The path of INTLAB toolbox and initialization.
% addpath("/path/to/your/INTLAB")
addpath("/Users/xfliu/App/Intlab_V14")

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

