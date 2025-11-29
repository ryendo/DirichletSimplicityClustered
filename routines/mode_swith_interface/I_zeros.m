function A = I_zeros(m, n)
    global INTERVAL_MODE;

    % Ensure dimensions are double (extract midpoint if intval)
    if isa(m, 'intval'), m = mid(m); end
    if isa(n, 'intval'), n = mid(n); end

    if INTERVAL_MODE
        A = intval(zeros(m, n));
    else
        A = zeros(m, n);
    end
end