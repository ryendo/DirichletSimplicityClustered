function val = f_tri_int(f, tri_intval)

    a = tri_intval(5); % 三角形の頂点 (a, b) の x 座標
    b = tri_intval(6); % 三角形の頂点 (a, b) の y 座標
    
    % 積分領域の設定
    xmin = 0;
    xmax = 1;
    ymin = @(x) (b/a) * (x - 1);
    ymax = @(x) (b/a) * x;
    
    % integral2 による数値積分
    val = integral2(f, xmin, xmax, ymin, ymax);

end
