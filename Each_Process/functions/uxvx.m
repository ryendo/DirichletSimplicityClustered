function output=uxvx(u,v,func_data)

    feature('setround', Inf);
    output_sup = sup(u)'*sup(func_data.Txx)*sup(v);
    feature('setround', -Inf);
    output_inf = inf(u)'*inf(func_data.Txx)*sup(v);
    feature('setround', 0.5);

    output = hull(output_inf,output_sup);
end