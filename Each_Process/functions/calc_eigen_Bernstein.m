function eigenValues = calc_eigen_Bernstein(af, bf)

    % Import Symbolic Math Toolbox
    syms x y a b

    % Degree of Bernstein polynomials
    n = 6;  % Adjust the degree here
    
    % Define the Bernstein basis functions
    L1 = (b*x - a*y);
    L2 = y;
    L3 = (1-a)*y + b*(x-a) - b*(1-a);
    bernstein_basis = sym(zeros(1, nchoosek(n, 2)));
    index = 1;
    for i = 1:n
        for j = 1:(n-i-1)
            k = n - i - j;
            bernstein_basis(index) = factorial(n) / (factorial(i) * factorial(j) * factorial(k)) * L1^i * L2^j * L3^k;
            index = index + 1;
        end
    end

    % Initialize mass matrix M and stiffness matrix K
    nBasis = length(bernstein_basis);
    M = sym(zeros(nBasis, nBasis));
    K = sym(zeros(nBasis, nBasis));

    % Compute mass matrix and stiffness matrix symbolically
    for i = 1:nBasis
        for j = 1:nBasis
            Ni = bernstein_basis(i);
            Nj = bernstein_basis(j);

            % Compute mass matrix
            % M(i,j) = ∫∫_Ω (N_i * N_j) dΩ
            M(i,j) = int(int(Ni * Nj, x, (a*y)/b, (-1+a)*(y-b)/b+a), y, 0, b);
            % Compute stiffness matrix
            % K(i,j) = ∫∫_Ω (∇N_i ⋅ ∇N_j) dΩ
            dNi_dx = diff(Ni, x);
            dNi_dy = diff(Ni, y);
            dNj_dx = diff(Nj, x);
            dNj_dy = diff(Nj, y);

            K(i,j) = int(int((dNi_dx * dNj_dx + dNi_dy * dNj_dy), x, (a*y)/b, (-1+a)*(y-b)/b+a), y, 0, b);
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

    % % Display sorted eigenvalues
    % disp('Sorted Eigenvalues:');
    % disp(eigenValues);
    % 
    % % Plot the first eigenfunction
    % % Create a grid for plotting
    % [X, Y] = meshgrid(linspace(0, 1, 100), linspace(0, bf, 100));
    % validRegion = inpolygon(X, Y, [0, 1, af], [0, 0, bf]);
    % 
    % % Calculate the first eigenfunction on the grid
    % eigenfunction = zeros(size(X));
    % for i = 1:nBasis
    %     eigenfunction = eigenfunction + eigenVectors(i, 2) * double(subs(subs(bernstein_basis(i), {a, b, x, y}, {af, bf, X, Y})));
    % end
    % 
    % % Set the invalid region to NaN
    % eigenfunction(~validRegion) = NaN;
    % 
    % % Plot the first eigenfunction
    % figure;
    % surf(X, Y, eigenfunction, 'EdgeColor', 'none');
    % colorbar;
    % title('First Eigenfunction');
    % xlabel('x');
    % ylabel('y');
    % zlabel('Eigenfunction value');
end
