function eigenValues = test_f_eigen_Bernstein(af, bf)
    llams = calc_eigen_bounds_cr(tri_intval,32);
    nu = llams(4);
    eig_m = 3;
    
    % Import Symbolic Math Toolbox
    syms x y
    a = tri_intval(5); b = tri_intval(6);

    mesh_int = get_mesh(tri_intval, 128);

    % Degree of Bernstein polynomials
    n = 5;  % Adjust the degree here
    bbasis = sym(zeros(1, nchoosek(n, 2)));
    nBasis = length(bbasis);

    % Check if files exist and load matrices
    % Define the Bernstein basis functions
    L1 = (b*x - a*y);
    L2 = y;
    L3 = (1-a)*y + b*(x-a) - b*(1-a);
    index = 1;
    for i = 1:n
        for j = 1:(n-i)
            k = n - i - j+1;
            bbasis(index) = factorial(n) / (factorial(i) * factorial(j) * factorial(k)) * L1^i * L2^j * L3^k;
            index = index + 1;
        end
    end


    % Initialize mass matrix M and stiffness matrix K
    M = intval(zeros(nBasis, nBasis));
    K = intval(zeros(nBasis, nBasis));

    % Compute mass matrix and stiffness matrix symbolically
    for i = 1:nBasis
        for j = i:nBasis  % Utilize symmetry by starting j from i
            [i, j];
            Ni = matlabFunction(bbasis(i)); % @(x,y)
            Nj = matlabFunction(bbasis(j));

            % Compute mass matrix
            % M(i,j) = ∫∫_Ω (N_i * N_j) dΩ
            M(i,j) = intval_int(Ni*Nj, tri_intval);

            % Compute stiffness matrix
            % K(i,j) = ∫∫_Ω (∇N_i ⋅ ∇N_j) dΩ
            
            dNi_dx = diff(Ni, x);
            dNi_dy = diff(Ni, y);
            dNj_dx = diff(Nj, x);
            dNj_dy = diff(Nj, y);

            K(i,j) = int(int((dNi_dx * dNj_dx + dNi_dy * dNj_dy), x, (a*y)/b, (-1+a)*(y-b)/b+a), y, 0, b);

            if i ~= j  % Fill the symmetric part
                M(j,i) = M(i,j);
                K(j,i) = K(i,j);
            end
        end
    end

    % Substitute numerical values of the vertices into symbolic matrices
    M_numeric = double(subs(subs(M, a, af), b, bf));
    K_numeric = double(subs(subs(K, a, af), b, bf));

    % Solve eigenvalue problem and sort eigenvalues
    [eigenVectors, eigenValues] = eig(K_numeric, M_numeric);
    eigenValues = diag(eigenValues);
    [eigenValues, sortIdx] = sort(eigenValues);
    eigenVectors = eigenVectors(:, sortIdx);

end
