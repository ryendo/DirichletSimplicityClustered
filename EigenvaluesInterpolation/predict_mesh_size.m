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
    lambda = v4(3)*0.995;

    % Required lower bound target:
    T = v1(2) + 0.2 * GAP;

    if T >= lambda
        h = 0.003;
        warning("Target lower bound T >= lambda; no positive h can satisfy the condition. h=0.002 is used.");
        return
    end

    % Solve formula:
    % h = (1/0.1893) * sqrt( lambda/T - 1 )
    coeff = 0.1893;
    h = (1/coeff) * sqrt(1/T - 1/lambda);

    % Make sure the lower bound of 3rd eigenvalue has less than 0.5%
    % relative error.
    % v(3) - v(3)/(1+(0.1893h)^2*v(3)) < v(3)*0.005
    % That is 1 + (0.1893h)^2*v(3) <  1/0.995, i.e., h <
    % sqrt((1-1/0.995)/v(3))/0.1893
    h2 = sqrt( (0.005/0.995) / v4(3) )/ 0.1893;

    if min(h2,h)<0.003
        fprintf("Predicted mesh size h by auto computation = %.10f\n", h,h2);
        
    end
    h = max(0.002,min(h2,h));

    fprintf("Predicted mesh size h = %.10f\n", h);
end
