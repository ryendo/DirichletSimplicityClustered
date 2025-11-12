% --- Minimal script: plot sup(plpu(1)) and sup(plpu(2)) only ---
% This keeps only what is necessary to compute SinvSinv^T, extract its two
% interval eigenvalues plpu(1:2), take their sups, and plot/save the curves.

% Requires INTLAB-style interval ops: I_intval, I_hull (or hull), I_veig, sup.

format long; tic

% Parameters
omega_N = 10;               % number of angular sectors on [0, pi/3]
ep      = I_intval('1E-5');   % perturbation radius
t       = I_hull(0, ep);      % interval parameter t in [0, ep]
I_PI    = I_intval('pi');     % interval version of pi

% Geometry of the reference equilateral triangle
x = I_intval('0.5');
y = sqrt(I_intval(3))/2;

% Output buffers for plotting
sup_plpu1 = zeros(omega_N,1);
sup_plpu2 = zeros(omega_N,1);
theta     = zeros(omega_N,1); % center angle (real) of each sector

% Create results folder if needed
if ~exist('results','dir'), mkdir results; end
out_png = ['results/plpu_sup_' datestr(now,'yyyy-mm-dd_HH-MM-SS') '2.png'];

% Main loop over sectors delta ∈ [ (i-1)/omega * pi/3 , i/omega * pi/3 ]
for idx = 1:omega_N
    idx
    % Interval angle for the current sector
    delta = I_hull((idx-1) * I_PI / (3*omega_N), idx * I_PI / (3*omega_N));
    delta = ((idx-1) * pi / (3*omega_N)+ idx * pi / (3*omega_N))/2

    % Direction (interval) d = (a, b) = (sin delta, -cos delta)
    a = sin(delta);
    b = -cos(delta);

    % Affine map components for the perturbed triangle
    tx = x + t * a;
    ty = y + t * b;

    % Inverse affine map S_inv and its symmetric product
    S_inv = [1, (x - tx) / ty; 0, y / ty];
    SinvSinvt = S_inv * S_inv'
    SinvSinvt = hull(SinvSinvt,SinvSinvt');
    % Make sure the enclosure is symmetric (INTLAB hull)
    % SinvSinvt = hull(SinvSinvt, SinvSinvt.')

    % Interval eigenvalues (two of them) of SinvSinvt
    plpu = I_veig(SinvSinvt, eye(2),1:2)

    % [V,D] = eig(I_mid(SinvSinvt))    
    % [mus,X] = verifyeig(SinvSinvt, D(1,1), V, eye(2,2))
    % % [mu3,X] = verifyeig(SinvSinvt, D(2,2), V(:,2), eye(2,2));
    % interval_mu = mus
    % [~,ind] = sort(I_mid(interval_mu));
    % plpu = interval_mu(ind);

    % Store the upper bounds (sups) for plotting
    sup_plpu1(idx) = inf(plpu(1))
    sup_plpu2(idx) = sup(plpu(2))

    % Real center angle of this sector for the x-axis
    theta(idx) = ((2*idx - 1) * pi) / (6 * omega_N);
end

% Plot
fig = figure('Color','w','Position',[200 200 900 500]);
plot(theta, sup_plpu1, '-+', 'LineWidth', 1.6); hold on;
plot(theta, sup_plpu2, '-+', 'LineWidth', 1.6);
grid on; box on;
xlabel('\delta (radian)');
ylabel('sup(plpu(i))');
legend({'sup(plpu(1))','sup(plpu(2))'}, 'Location', 'best');
title('sup(plpu(1)) and sup(plpu(2)) vs. direction \delta');

% Save PNG
try
    exportgraphics(fig, out_png, 'Resolution', 300);
catch
    saveas(fig, out_png); % fallback for older MATLAB
end

toc
