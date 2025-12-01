function h = predict_mesh_size(region_cell)
%   let v1(:) be the eigenvalues at (a1,b1)
%   let v4(:) be the eigenvalues at (a4,b4)
%   GAP = v1(3)-v1(2)
%   lambda = v4(3)
%   LB = lambda/(1+(0.1893*h)^2*lambda)
%   Require LB - v1(2) >= 0.2 * GAP
%   Solve for h.

    x1 = region_cell(1);
    x2 = region_cell(2);
    t1 = region_cell(3);
    t2 = region_cell(4);
    a1 = x1;  b1 = x1*tan(t1);
    a4 = x2;  b4 = x2*tan(t2);

    % --- Compute eigenvalues at both points ---
    [e1,e2,e3,e4] = get_approximate_eigenvalue(a1,b1);
    v1 = [e1, e2, e3, e4];

    [e1,e2,e3,e4] = get_approximate_eigenvalue(a4,b4);
    v4 = [e1, e2, e3, e4];

    % --- Data ---
    GAP = v1(3) - v1(2);
    lambda = v4(3);

    % Required lower bound target:
    T = v1(2) + 0.18 * GAP;

    if T >= lambda
        h = 0.002;
        warning("Target lower bound T >= lambda; no positive h can satisfy the condition. h=0.002 is used.");
        return
    end

    % Solve formula:
    % h = (1/0.1893) * sqrt( lambda/T - 1 )
    coeff = 0.1893;
    h = (1/coeff) * sqrt(1/T - 1/lambda);

    fprintf("Predicted mesh size h = %.10f\n", h);
end