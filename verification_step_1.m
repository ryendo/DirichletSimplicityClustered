function interval_mu = verification_step_1(t, delta, mat, u1, u2, u3)
% RUN_ALGO1_STEP Computes rigorous bounds for difference quotients.
%   interval_mu = verification_step_1(t, delta, mat, u1, u2, u3)
%
%   Inputs:
%       t     : Interval [0, ep] representing perturbation magnitude.
%       delta : Interval representing the perturbation angle sector.
%       mat   : Struct containing FEM matrices (xx, xy, yy, grad, l2).
%       u1,u2,u3 : Basis vectors (approximate eigenfunctions).
%
%   Output:
%       interval_mu : [mu2, mu3] rigorous interval bounds.
%
%   This function encapsulates the logic of "Algorithm 1" from the paper.

    % --- Core Algorithm 1 Logic (unchanged) ---
    
    % Calculate the perturbation direction 'd' and the matrix 'Pt' (P_t^e).
    d = d_a(delta); a = d(1); b = d(2);
    x = I_intval('0.5'); tx = x + t * a; 
    y = sqrt(I_intval('3')) / 2; ty = y + t * b;
    Pt = [t * a / (ty^2), -a * y / (ty^2); -a * y / (ty^2), -b * (y + ty) / (ty^2)];

    % [LIU] tx: \tilde{x}, ty: \tilde{y} 
    S = [1, (tx-x) / y; 0, ty / y];
    % Removed the square from norm_StStT related to Lemma 3.7
    norm_StStT = I_intval(I_sup(norm(S*S',2)));
    
    S_inv = [1, (x-tx) / ty; 0, y / ty];
    SinvSinvt = S_inv * S_inv';
    
    % Removed the square from norm_StStT_inv related to Lemma 3.7
    norm_StStT_inv = I_intval(I_sup(norm(SinvSinvt,2))); 
    
    % Define the exact, known eigenvalues for the equilateral triangle.
    lam = [I_intval('16') / I_intval('3') * I_pi^2, ...
           I_intval('112') / I_intval('9') * I_pi^2, ...
           I_intval('112') / I_intval('9') * I_pi^2, ...
           I_intval('64') / I_intval('3') * I_pi^2];
    lam1 = lam(1); lam2 = lam(2); lam3 = lam(3); lam4 = lam(4);
    
    rho = lam(4);
    
    % Calculate the Rayleigh quotients for the approximate eigenfunctions.
    lam1h = I_intval(I_sup(grad(mat, u1, u1) / L2(mat, u1, u1))); 

    matA_lam2hlam3h = [grad(mat, u2, u2), grad(mat, u2, u3); grad(mat, u3, u2), grad(mat, u3, u3)];
    matB_lam2hlam3h = [L2(mat, u2, u2), L2(mat, u2, u3); L2(mat, u3, u2), L2(mat, u3, u3)];
    lam2hlam3h = norm(matB_lam2hlam3h \ matA_lam2hlam3h,2);
    lam3h = max(lam2hlam3h);

    % Upper bound of \tilde{\lambda}_N using the estimation in Lemma 3.7.
    lam1_tilde =  lam1*(norm_StStT_inv * norm_StStT);
    lam3_tilde =  lam3*(norm_StStT_inv * norm_StStT);
    
    % --- Rigorous Error Analysis Implementation ---
    matF = [L2(mat, u1, u2); L2(mat, u1, u3)];
    matG = L2(mat, u1, u1);
    matH = [L2(mat, u2, u2), L2(mat, u2, u3); L2(mat, u3, u2), L2(mat, u3, u3)];
    etaF = norm(matF * matF', 2);
    etaG = norm(I_eye(1, 1) - matG, 2);
    etaH = norm(I_eye(2, 2) - matH, 2);
    epb_hat12 = etaF / ((1 - etaG) * (1 - etaH));
    
    dbE1E1h = sqrt((lam1h - lam1) / (rho - lam1h));
    thb2 = (rho - lam2) * (epb_hat12 + dbE1E1h)^2;
    dbE2E2h = sqrt((lam3h - lam2 + thb2) / (rho - lam2));
    dbbarE2E2h = sqrt(2 - 2 * sqrt(1 - dbE2E2h^2));
    
    dbE1E1t = sqrt((lam1_tilde - lam1) / (rho - lam1));
    thb2t = (rho - lam2) * dbE1E1t^2;
    dbE2E2t = sqrt((lam3_tilde - lam2 + thb2t) / (rho - lam2));
    dbbarE2E2t = sqrt(2 - 2 * sqrt(1 - dbE2E2t^2));
    
    % Calculate the error terms 'eta_hat' and 'eta_tilde' from Lemma 3.4
    eta_hat = sqrt(lam3 * dbbarE2E2h^2 + lam3h - lam2);
    eta_tilde = sqrt(lam3 * dbbarE2E2t^2 + lam3_tilde - lam2);

    % Calculate the final error bounds 'Err_F' and 'Err_b'
    Err_F = I_intval(I_sup(norm(Pt, 2) * (2 * sqrt(lam3h) * eta_hat + sqrt(lam3) * eta_tilde)))
    Err_b = I_intval(I_sup(2 * dbbarE2E2h + dbbarE2E2t))

    % Matrix M
    M = [F(mat, Pt, u2, u2), F(mat, Pt, u2, u3); F(mat, Pt, u3, u2), F(mat, Pt, u3, u3)];
    N = [L2(mat, u2, u2), L2(mat, u2, u3); L2(mat, u3, u2), L2(mat, u3, u3)];

    % --- Construct Interval Matrices ---
    M_ = M + I_hull(-Err_F, Err_F);
    N_ = N + I_hull(-Err_b, Err_b);
    
    % --- Solve the Generalized Eigenvalue Problem ---
    [V,D] = eig(I_mid(M_), I_mid(N_));
    [mu2_,~] = verifyeig(M_, D(1,1), V(:,1), N_);
    [mu3_,~] = verifyeig(M_, D(2,2), V(:,2), N_);

    mus =  I_eigs(M_, N_, 2, 'sm');
    
    interval_mu = [mu2_, mu3_];
    [~,ind] = sort(I_mid(interval_mu));
    interval_mu = interval_mu(ind);

end

% --- Helper Functions (Local) ---

function d = d_a(delta)
    d = [sin(delta), -cos(delta)];
end

function f = F(mat, Pt, u, v)
    uxvx = u' * mat.xx * v; uyvy = u' * mat.yy * v;
    uxvy = u' * mat.xy * v; uyvx = v' * mat.xy * u;
    f = Pt(1, 1) * uxvx + Pt(1, 2) * uyvx + Pt(2, 1) * uxvy + Pt(2, 2) * uyvy;
end

function f = L2(mat, u, v)
    f = u' * mat.l2 * v;
end

function f = grad(mat, u, v)
    f = u' * mat.grad * v;
end