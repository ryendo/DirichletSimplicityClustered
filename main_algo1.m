% --- Main Script for Step 1: Rigorous Difference Quotient Estimation ---
% This script implements the computer-assisted proof for the Omega_up region,
% i.e., for triangles very close to equilateral. It computes guaranteed bounds
% for the difference quotients of the 2nd and 3rd eigenvalues to prove they
% are separate. This corresponds to Algorithm 1 in the paper.

% Configure the INTLAB library for interval arithmetic and set display format.
% All subsequent calculations involving intervals will have rigorous error bounds.
% my_intlab_mode_config;
tic
format compact short infsup

% --- Set main parameters for the computation ---
% omega_N: Number of angular intervals to discretize the perturbation direction.
% mesh_N:  Parameter controlling the fineness of the finite element mesh.
% ep:      Maximum perturbation distance 't' from the equilateral triangle.
% ord:     The polynomial degree (order) of the finite element method.
% t:       The interval [0, ep] representing the perturbation magnitude.

omega_N = 1000; % 1000
mesh_N = 8; % 8
ep = I_intval('1E-5'); % (1E-5): (ord=5,mesh_N=8,omega_N = 1000)
ord = 5; % 5: Order of FEM
t = I_hull(0,ep); % Interval parameter

% Define the base (equilateral) triangle using interval arithmetic.
tri_I_intval = [I_intval('0'), I_intval('0'); I_intval('1'), I_intval('0'); 1/I_intval(2), sqrt(I_intval(3))/2];

% --- Initial Finite Element Method (FEM) Calculation ---
% Perform the FEM analysis on the equilateral triangle to obtain:
% - lamhs: Approximate eigenvalues.
% - inner_uhat: Approximate eigenfunctions (as coefficient vectors).
% - A_...: The discretized system matrices (stiffness, mass, etc.).
[lamhs, inner_uhat, A_inner_grad, A_inner_L2, A_inner_xx, A_inner_xy, A_inner_yy, bd_dof_idx] = calc_eigen_bounds_any_order_for_quotients(tri_I_intval, mesh_N, ord);

% Bundle the FEM matrices into a struct for convenient access.
mat = struct('xx', A_inner_xx, 'xy', A_inner_xy, 'yy', A_inner_yy, 'grad', A_inner_grad, 'l2', A_inner_L2);

% --- Orthonormalize the basis for the approximate eigenspace ---
% The theory requires an orthonormal basis for the eigenspace corresponding to
% the multiple eigenvalue (lambda_2 = lambda_3).
% u1, u2, u3 are the approximate eigenfunctions.
% The Gram-Schmidt process is used to ensure u3 is orthogonal to u2.
u1 = inner_uhat(:, 1); 
u2 = inner_uhat(:, 2); 
u3 = inner_uhat(:, 3);
u1 = u1 / sqrt(L2(mat, u1, u1));
u2 = u2 / sqrt(L2(mat, u2, u2));
u3 = u3 - u2 * L2(mat, u3, u2); 
u3 = u3 / sqrt(L2(mat, u3, u3));

file_name = ['results/quotients_' datestr(now, 'yyyy-mm-dd_HH-MM-SS') '.csv'];

% --- Main loop to compute difference quotients for all directions ---
% Iterates through 'omega_N' angular sectors to cover all perturbation
% directions from the equilateral triangle.
% for idx = 1:omega_N
for idx = [1,omega_N]
    % Define the interval 'delta' for the current perturbation angle sector.
    delta = I_hull((idx-1) * I_pi / (3 * omega_N), idx * I_pi / (3 * omega_N));

    % Call the core function to calculate the difference quotient bounds.
    [approx_mu, interval_mu] = quotient_calc(t, delta, mat, u1, u2, u3); % Compute difference quotients

    % Store the results: index and the rigorous interval bounds for the
    % difference quotients of the 2nd and 3rd eigenvalues.
    info_mat = [idx, I_inf(interval_mu(1)), I_sup(interval_mu(1)), I_inf(interval_mu(2)), I_sup(interval_mu(2))];

    % Append the newly computed bounds to the CSV file.
    writematrix(info_mat, file_name, 'WriteMode', 'append'); % Write results to file
end

% --- Helper function to calculate the 2D direction vector 'e'. ---
% This corresponds to the vector e(delta) defined in Figure 6 of the paper.
function d = d_a(delta)
    d = [sin(delta), -cos(delta)];
end

% --- Helper function to compute the bilinear form F_P(u, v). ---
% Calculates (P * grad(u), grad(v)) using the pre-computed FEM matrices.
% This form is central to the difference quotient formula (Theorem 3.2).
function f = F(mat, Pt, u, v)
    uxvx = u' * mat.xx * v; uyvy = u' * mat.yy * v;
    uxvy = u' * mat.xy * v; uyvx = v' * mat.xy * u;
    f = Pt(1, 1) * uxvx + Pt(1, 2) * uyvx + Pt(2, 1) * uxvy + Pt(2, 2) * uyvy;
end

% --- Helper function to compute the L2 inner product (u, v). ---
function f = L2(mat, u, v)
    f = u' * mat.l2 * v;
end

% --- Helper function to compute the gradient inner product (grad(u), grad(v)). ---
function f = grad(mat, u, v)
    f = u' * mat.grad * v;
end

% --- Core function to compute guaranteed bounds on eigenvalue difference quotients. ---
% This function implements the main steps of Algorithm 1 from the paper.
function [approx_mu, interval_mu] = quotient_calc(t, delta, mat, u1, u2, u3)

    % Calculate the perturbation direction 'd' and the matrix 'Pt' (P_t^e in the paper).
    d = d_a(delta); a = d(1); b = d(2);
    x = I_intval('0.5'); tx = x + t * a; 
    y = sqrt(I_intval(3)) / 2; ty = y + t * b;
    Pt = [t * a / (ty^2), -a * y / (ty^2); -a * y / (ty^2), -b * (y + ty) / (ty^2)];

    % Estimate eigenvalue bounds on the perturbed domain ('lam_tilde')
    % by calculating perturbation factors from the affine map 'S'.
    % S = [1, (tx-x) / y; 0, ty / y];
    % plpu_ = I_veig(hull(S*S',(S*S')'), eye(2, 2), 1:2);
    % pu_ = I_intval(I_sup(plpu_(2)));
    % S_inv = [1, (x-tx) / ty; 0, y / ty];
    % SinvSinvt = S_inv * S_inv'; SinvSinvt=hull(SinvSinvt,SinvSinvt');
    % plpu = I_veig(SinvSinvt, eye(2, 2), 1:2);
    % pl = I_intval(I_inf(plpu(1))); pu = I_intval(I_sup(plpu(2)));

    % [LIU] tx: \tilde{x}, ty: \tilde{y} 
    % [LIU] pl, pu: lower and upper bounds for the largest eigenvalue of S*S'.
    S = [1, (tx-x) / y; 0, ty / y];
    pu_ = I_intval(I_sup(norm(S*S',2)^2));   %[LIU Request] Are both pu_ and pu are needed? 
    S_inv = [1, (x-tx) / ty; 0, y / ty];
    SinvSinvt = S_inv * S_inv';
    pl = I_intval(I_inf(1/norm(S*S',2)^2));
    pu = I_intval(I_sup(norm(SinvSinvt,2)^2)); 
    
    % Define the exact, known eigenvalues for the equilateral triangle.
    % These serve as the base values for the difference quotient calculation.
    lam = [I_intval('16') / I_intval('3') * I_intval('pi')^2, I_intval('112') / I_intval('9') * I_intval('pi')^2, I_intval('112') / I_intval('9') * I_intval('pi')^2, I_intval('64') / I_intval('3') * I_intval('pi')^2];
    lam1 = lam(1); lam2 = lam(2); lam3 = lam(3); lam4 = lam(4);
    
    % [LIU Request] "the shift parameter 'rho'" -> "the lower eigenvalue bound 'rho'" 
    % Define the shift parameter 'rho' for the error estimation formulas.
    % It must satisfy lambda_N < rho <= lambda_{N+1}. Here, N=3.
    rho = lam(4) * pl;
    
    % Calculate the Rayleigh quotients for the approximate eigenfunctions.
    % These are the approximate eigenvalues (e.g., lambda_{N,h} in the theory).
    % [LIU] Only rigorous upper bound of lambda_h is needed.
    lam1h = I_intval(I_sup(grad(mat, u1, u1) / L2(mat, u1, u1))); 
    lam3h = I_intval(I_sup(grad(mat, u3, u3) / L2(mat, u3, u3)));

    %  [LIU] Upper bound of \tilde{\lambda}_N using the estimation in Lemma 3.7.
    %  [LIU Request] Upon Lemma 3.7, the 2-norm of matrix Sinv*Sinvt is used. But the definition of pu_ takes the square of the matrix.
    lam1_tilde =  (lam1*pu) * pu_;
    lam3_tilde = (lam3*pu) * pu_;
    
    % --- Rigorous Error Analysis Implementation ---
    % This block directly implements the error estimation theory from Section 3.3
    % and Lemma 2.1 to bound the errors between different eigenspaces.
    % [LIU] The estimation of \epsion_{a}^h, \epsilon_{b}^h refers to Lemma 2 of [11].
    % dbbarE2E2h: Error between the true eigenspace E and the computed approximate eigenspace \hat{E}.
    % dbbarE2E2t: Error between the true eigenspace E  and the mapped perturbed eigenspace \tilde{E}_t.
    matF = [L2(mat, u1, u2); L2(mat, u1, u3)];
    matG = L2(mat, u1, u1);
    matH = [L2(mat, u2, u2), L2(mat, u2, u3); L2(mat, u3, u2), L2(mat, u3, u3)];
    etaF = norm(matF * matF', 2);
    etaG = norm(eye(1, 1) - matG, 2);
    etaH = norm(eye(2, 2) - matH, 2);
    epb_hat12 = etaF / ((1 - etaG) * (1 - etaH));
    
    dbE1E1h = sqrt((lam1h - lam1) / (lam4 - lam1h));
    thb2 = (lam4 - lam2) * (epb_hat12 + dbE1E1h)^2;
    dbE2E2h = sqrt((lam3h - lam2 + thb2) / (lam4 - lam2));
    dbbarE2E2h = sqrt(2 - 2 * sqrt(1 - dbE2E2h^2));
    
    dbE1E1t = sqrt((lam1_tilde - lam1) / (rho - lam1));
    thb2t = (rho - lam2) * dbE1E1t^2;
    dbE2E2t = sqrt((lam3_tilde - lam2 + thb2t) / (rho - lam2));
    dbbarE2E2t = sqrt(2 - 2 * sqrt(1 - dbE2E2t^2));
    
    % Calculate the error terms 'eta_hat' and 'eta_tilde' from Lemma 3.4,
    % which bound the error in the gradient norm.
    eta_hat = sqrt(lam3 * dbbarE2E2h^2 + lam3h - lam2);
    eta_tilde = sqrt(lam3 * dbbarE2E2t^2 + lam3_tilde - lam2);

    % Calculate the final error bounds 'Err_F' and 'Err_b' for the
    % system matrices, as defined in Lemma 3.4.
    Err_F = I_intval(I_sup(norm(Pt, 2) * (2 * sqrt(lam3h) * eta_hat + sqrt(lam3) * eta_tilde)))
    Err_b = I_intval(I_sup(2 * dbbarE2E2h + dbbarE2E2t))

    % Compute the approximate matrices 'M' (hat{M}_t) and 'N' (hat{N}_t)
    % for the generalized eigenvalue problem of the difference quotient.
    % M = [F(mat, Pt, u2, u2), F(mat, Pt, u2, u3); F(mat, Pt, u3, u2), F(mat, Pt, u3, u3)]; 
    M = [F(mat, Pt, u3, u3), F(mat, Pt, u3, u2); F(mat, Pt, u2, u3), F(mat, Pt, u2, u2)]
    N = [L2(mat, u2, u2), L2(mat, u2, u3); L2(mat, u3, u2), L2(mat, u3, u3)]

    % --- Construct Interval Matrices ---
    % This is the key step for rigorous computation. It adds the guaranteed
    % error bounds (Err_F, Err_b) to the approximate matrices M and N
    % to create interval matrices M_ and N_ that are guaranteed to contain the true ones.
    % M = I_hull(M,M'); N = I_hull(N,N');

    M_ = M + I_hull(-Err_F, Err_F);
    N_ = N + I_hull(-Err_b, Err_b);
    
    % % --- Solve the Generalized Eigenvalue Problem ---
    % % 'approx_mu': The eigenvalues of the approximate problem.
    % % 'interval_mu': The eigenvalues of the interval problem, which gives the
    % % final, guaranteed bounds on the difference quotients.    

    % approx_mu = I_veig(M, N, 1:2)
    % interval_mu = I_veig(M_, N_, 1:2)

    % --- Solve the Generalized Eigenvalue Problem ---

    [V,D] = eig(I_mid(M),I_mid(N));    
    [mu2,X] = verifyeig(M, D(1,1), V(:,1),N);
    [mu3,X] = verifyeig(M, D(2,2), V(:,2),N);
    approx_mu = [mu2,mu3];
    [~,ind] = sort(I_mid(approx_mu));
    approx_mu = approx_mu(ind)

    [V,D] = eig(I_mid(M_),I_mid(N_));    
    [mu2_,X] = verifyeig(M_, D(1,1), V(:,1),N_);
    [mu3_,X] = verifyeig(M_, D(2,2), V(:,2),N_);
    interval_mu = [mu2_,mu3_]
    [~,ind] = sort(I_mid(interval_mu));
    interval_mu = interval_mu(ind)

    toc
end
