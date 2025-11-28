% Verified eigenvalue enclosures via similarity + Gershgorin.
% Assumes that V from eig(mid(A)) is nonsingular. Then X = V\(A*V) ⊇ V^{-1}AV.
% Row Gershgorin: r upper-bounds off-diagonal row sums for all X*∈X ⇒
% eig_vals = diag(X) ± r encloses the spectrum (disjoint boxes certify counts).


function eig_vals = veig_ge(A)

    A = I_intval(A);
    [V D] = eig(I_mid(A));
    X = I_intval(V) \ (A*I_intval(V));
    r = sum(abs(X), 2) - abs(diag(X));
    eig_vals = diag(X)+I_hull(-r,r);
    [~,id] = sort(I_mid(eig_vals));
    eig_vals = eig_vals(id);

    eig_vals_before = eig_vals(id)

    %----------------------
    A = I_intval(A);
    [V D] = eig(I_mid(A));
    X = I_intval(V) \ (A*I_intval(V));

    n = size(X,1);
    r = I_zeros(n,1);    
    for i = 1:n
        for j = [1:i-1, i+1:n]
            r(i) = r(i) + abs(X(i,j));
        end
    end
    eig_vals = diag(X)+I_hull(-r,r);
    [~,id] = sort(I_mid(eig_vals));
    eig_vals = eig_vals(id);

    eig_vals_after = eig_vals(id)

    %----------------------

    [V,D] = eig(I_mid(A));    
    [mu2,X] = verifyeig(A, D(1,1), V(:,1));
    [mu3,X] = verifyeig(A, D(2,2), V(:,2));
    eig_vals_after2 = [mu2,mu3]

end