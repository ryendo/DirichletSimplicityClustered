%clear;clc;
mesh_path   = './mesh_data/UnitSquare8x8/';
%output_path = 'C:\Users\youch\Documents\MATLAB\SteklovEigen\lower_bound_gather\output\';
vert = load([mesh_path, 'vert.dat']);
edge = load([mesh_path, 'edge.dat']);
tri  = load([mesh_path, 'tri.dat']);
bd   = load([mesh_path, 'bd.dat']);

%write_mesh_xml(vert, tri, [mesh_path, 'usd.xml'])

% trimesh(tri, vert(:,1), vert(:,2), 'color', 'k')
% axis equal
% axis off

neig = 4;
mesh_size = [size(vert, 1), size(tri, 1)]

% %================================================================
% disp('compute steklov eig using CR element');
% [CR_eig, Ch] = CR_steklov_eig(vert, edge, tri, bd, neig);
% Ch = Ch(1);
% disp('compute steklov eig lower bound by our thm')
% CR_eig_low = CR_eig ./ (1 + CR_eig .* Ch^2);

%================================================================
disp('compute steklov eig using Lagrange element');
Lagrange_order = 2;
[eig_value, eig_func, ~, ~] = Lagrange_dirichlet_eig(Lagrange_order, vert, edge, tri, bd, neig);


% %================================================================
% disp('compute steklov eig lower bound by L-G thm')
% A0 = LA_eigf' * LA_A * LA_eigf;
% A1 = LA_eigf' * LA_M * LA_eigf;
% rho = 2.5;
% 
% RT_order = Lagrange_order;
% RT_femf = [];
% for i = 1:neig
%     f = f_for_Hdiv_problem{i, 1};
%     [femf, ~, ~, RT_A] = RT_Hdiv_problem(RT_order, vert, edge, tri, bd, f);
%     RT_femf = [RT_femf, femf];
% end
% A_lg = RT_femf' * RT_A * RT_femf;
% 
% %A2_diff = A_lg - A2
% 
% AL = A0 - rho * A1;
% BL = A0 - 2*rho*A1 + rho*rho*A_lg;
% 
% [~, LG_eig_low] = eig(AL, BL);
% LG_eig_low = sort(diag(LG_eig_low));
% LG_eig_low = rho - rho./(1-LG_eig_low(end:-1:1));
% 
% eig_cr_crl_la_lgl = [CR_eig, CR_eig_low, LA_eig, LG_eig_low]
