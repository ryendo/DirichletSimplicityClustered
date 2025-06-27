format compact short infsup

a = intval('0.5');
b = intval('5')/1000;

tri_intval=[I_intval('0'), I_intval('0'); I_intval('1'), I_intval('0'); a,b];

isLG = 1;
ord = 1;
N_rho = 64;
N_LG = 8;

lam1=calc_eigen_bounds_any_order_lam1(tri_intval,N_LG,N_rho,ord,isLG);

peri_tri = 1+sqrt(a^2+b^2)+sqrt((a-1)^2+b^2);
area_tri = b/2;

F = lam1*area_tri-I_pi^2/16*peri_tri^2/area_tri
F_opt = 7*sqrt(intval(3))*I_pi^2/12