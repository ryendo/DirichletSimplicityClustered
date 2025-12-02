function draw_regions(regions)
    region_n = size(regions,1);
    for k = 1:region_n
          x1 = regions(k,1);
          x2 = regions(k,2);
          t1 = regions(k,3);
          t2 = regions(k,4);
          a1 = x1;  b1 = x1*tan(t1);
          a2 = x1;  b2 = x1*tan(t2);
          a3 = x2;  b3 = x2*tan(t1);
          a4 = x2;  b4 = x2*tan(t2);
          hold on
          plot([a1,a3,a4,a2,a1],[b1,b3,b4,b2,b1],'-');
    end