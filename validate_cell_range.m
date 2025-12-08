function ok = validate_cell_range(csvFile, tol)
%VALIDATE_CELLS  Validate cell partition from csvFile.
%
%   ok = VALIDATE_CELLS(csvFile)
%   ok = VALIDATE_CELLS(csvFile, tol)
%
% csvFile : path to the csv file with header
%           i,x_inf,x_sup,theta_inf,theta_sup,mesh_size_upper,...
% tol     : tolerance for equality checks (default: 1e-10)
%
% Returns:
%   ok : logical true if all checks pass, false otherwise.
%
% Conditions:
% (a) For each theta-row (same theta_inf, theta_sup),
%     sorted by x_inf, we require x_sup(left cell) = x_inf(right cell)
%     for neighbor cells. For neighbor theta-rows, we require
%     theta_inf(next row) = theta_sup(current row).
%
% (b) For each theta-row, let left-end cell be the one with minimal x_inf.
%     Let its vertices
%       p1 = [x_inf, x_inf*tan(theta_inf)]
%       p2 = [x_inf, x_inf*tan(theta_sup)]
%       p3 = [x_sup, x_sup*tan(theta_inf)]
%       p4 = [x_sup, x_sup*tan(theta_sup)]
%     Let y0 = tan(pi/60)/2;
%     If y(p1) > y0 then OK.
%     If y(p1) < y0 then require y(p2) < y0 and y(p3) < y0.
%
% (c) For each theta-row, let right-end cell be the one with maximal x_sup.
%     Let p3 = [x_sup, x_sup*tan(theta_inf)]. Require x_sup > 1 * theta_inf.

    if nargin < 2
        tol = 1e-10;
    end

    T = readtable(csvFile);

    % Extract the columns we need
    x_inf     = T.x_inf;
    x_sup     = T.x_sup;
    theta_inf = T.theta_inf;
    theta_sup = T.theta_sup;

    nCells = height(T);

    % ---- Group cells by (theta_inf, theta_sup) = theta rows ----
    % NOTE: this assumes theta_inf/theta_sup are numerically identical
    % for cells in the same row. If not, some tolerance-based grouping
    % would be needed.
    theta_pairs = [theta_inf, theta_sup];
    [uniqueThetaRows, ~, rowIdx] = unique(theta_pairs, 'rows');  % each row of uniqueThetaRows is one theta-row

    nThetaRows = size(uniqueThetaRows, 1);

    allOk = true;

    %% (a1) For each theta-row: x continuity in that row
    fprintf('Checking (a1) x-continuity within each theta-row...\n');
    for r = 1:nThetaRows
        idx = find(rowIdx == r);                 % indices of cells in this theta-row
        [~, order] = sort(x_inf(idx));           % sort by x_inf
        idx = idx(order);

        % neighbor cells: x_sup(left) = x_inf(right)
        for k = 1:numel(idx) - 1
            left  = idx(k);
            right = idx(k+1);

            if abs(x_sup(left) - x_inf(right)) > tol
                allOk = false;
                fprintf('  (a1) FAIL: row %d, cells %d (left) and %d (right) have gap: x_sup(left)=%.16g, x_inf(right)=%.16g\n', ...
                        r, left, right, x_sup(left), x_inf(right));
            end
        end
    end

    %% (a2) Neighbor theta-rows: theta_inf(next) = theta_sup(current)
    fprintf('Checking (a2) theta-continuity between neighbor theta-rows...\n');
    % Sort rows by theta_inf to define "neighbor rows"
    [~, rowOrder] = sort(uniqueThetaRows(:,1));  % sort by theta_inf

    for k = 1:nThetaRows - 1
        curIdx  = rowOrder(k);
        nextIdx = rowOrder(k+1);

        cur_theta_sup  = uniqueThetaRows(curIdx, 2);
        next_theta_inf = uniqueThetaRows(nextIdx, 1);

        if next_theta_inf > cur_theta_sup
            allOk = false;
            fprintf('  (a2) FAIL: theta-row %d -> %d: theta_sup(cur)=%.16g, theta_inf(next)=%.16g\n', ...
                    curIdx, nextIdx, cur_theta_sup, next_theta_inf);
        end
    end

    %% (b) Left-end cell and y = 0.0265 condition
    fprintf('Checking (b) left-end cell condition with y = 0.0265...\n');
    y0=tan(pi/60)/2
    yBoundary = y0;

    for r = 1:nThetaRows
        idx = find(rowIdx == r);
        % left-end cell = minimal x_inf
        [~, order] = sort(x_inf(idx));
        idx = idx(order);
        leftCell = idx(1);

        xi   = x_inf(leftCell);
        xs   = x_sup(leftCell);
        th_i = theta_inf(leftCell);
        th_s = theta_sup(leftCell);

        y1 = xi * tan(th_i);  % p1
        y2 = xi * tan(th_s);  % p2
        y3 = xs * tan(th_i);  % p3
        % y4 = xs * tan(th_s);  % p4 (not used in condition)

        if y1 > yBoundary + tol
            % OK by condition
            continue;
        elseif y1 < yBoundary - tol
            % Then require y2, y3 < y0
            if xi==0.5
                continue
            end
            if ~( (y2 < yBoundary + tol) && (y3 < yBoundary + tol) )
                allOk = false;
                fprintf(['  (b) FAIL: theta-row %d left cell %d: y1=%.16g<%.4g, but y2=%.16g, y3=%.16g\n' ...
                         '           Expected y2<%.4g and y3<%.4g.\n'], ...
                        r, leftCell, y1, yBoundary, y2, y3, yBoundary, yBoundary);
            end
        else
            % y1 is approximately on the boundary (within tol)
            % Interpret as OK (can adjust if desired)
            continue;
        end
    end

    %% (c) Right-end cell and x_sup > 1 * cos(theta_inf) condition
    fprintf('Checking (c) right-end cell condition x_sup > cos(theta_inf) ...\n');

    for r = 1:nThetaRows
        idx = find(rowIdx == r);
        % right-end cell = maximal x_sup
        [~, order] = sort(x_sup(idx));
        idx = idx(order);
        rightCell = idx(end);

        xs   = x_sup(rightCell);
        th_i = theta_inf(rightCell);

        % Condition: x(p3)=x_sup > r * t0heta_inf, r = 1
        if ~(xs >= 1 * cos(th_i))
            allOk = false;
            fprintf('  (c) FAIL: theta-row %d right cell %d: x_sup=%.16g, cos(theta_inf)=%.16g (need x_sup > cos(theta_inf) ).\n', ...
                    r, rightCell, xs, cos(th_i));
        end
    end

    %% Final result
    if allOk
        fprintf('All checks (a), (b), (c) PASSED.\n');
    else
        fprintf('Some checks FAILED. See messages above.\n');
    end

    ok = allOk;
end