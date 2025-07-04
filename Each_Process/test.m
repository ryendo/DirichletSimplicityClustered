

upper_right = [I_intval('0'), I_intval('0'); I_intval('1'), I_intval('0'); intval('0.5'), sqrt(intval(3))/2-10^(-5)];
N_LG = 4; % The Lehmann-Goerisch method is set to compute bounds for the first 4 eigenvalues.
fem_ord = 5; % A 5th-degree finite element method is used for high accuracy.
llams = calc_eigen_bounds_any_order(upper_right, N_LG, 64, fem_ord, 1)
                    