

tri_intval = [I_intval('0'), I_intval('0'); I_intval('1'), I_intval('0'); 1/intval('2'), sqrt(intval('3'))/2-10^(-5)];
N_LG = 4;
N_rho = 64;
fem_ord = 5;
isLG = 1;

% tri_intval = [I_intval('0'), I_intval('0'); I_intval('1'), I_intval('0'); 1/intval('2'), sqrt(intval('3'))/2];
% N_LG = 4;
% N_rho = 32;
% fem_ord = 5;
% isLG = 0;


calc_eigen_bounds_any_order(tri_intval,N_LG,N_rho,fem_ord,isLG)