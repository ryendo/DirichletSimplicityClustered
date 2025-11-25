%The path of the codes for switch between verified computing and approximate computing.
addpath('./lib_eigenvalue_bound/');
addpath('./lib_mesh/');
addpath('./mode_switch_interface/');  
addpath("./veigs/");

%The path of INTLAB toolbox and initialization.
addpath('~/app/Intlab_V12/');

clear INTERVAL_MODE;
global INTERVAL_MODE;

%INTERVAL_MODE=1; for rigorous computing based on interval arithmetic.
%INTERVAL_MODE=0; for approximate computing with rounding error inside.
INTERVAL_MODE=1;

if INTERVAL_MODE == 1
    startintlab;
end