function nodes=get_xy(region_cell)
    x1 = region_cell(1);
    x2 = region_cell(2);
    t1 = region_cell(3);
    t2 = region_cell(4);
    a1 = x1;  b1 = x1*tan(t1);
    a2 = x1;  b2 = x1*tan(t2);
    a3 = x2;  b3 = x2*tan(t1);
    a4 = x2;  b4 = x2*tan(t2);
    nodes=[a1,b1;a2,b2;a3,b3;a4,b4];
end