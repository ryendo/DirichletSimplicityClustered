
function val = Fr_t(u,v,p,t)
    
    t_p = p+t*[-1,0];
    r = p(1); theta = p(2);
    t_r = t_p(1); t_theta= t_p(2);

    S = [1,(t_r*cos(t_theta)-r*cos(theta))/(r*sin(theta));0,t_r*sin(t_theta)/(r*sin(theta))];
    Pt = (S^(-1)*transpose(S^(-1))-intval(eye(2,2)))/t

    

end