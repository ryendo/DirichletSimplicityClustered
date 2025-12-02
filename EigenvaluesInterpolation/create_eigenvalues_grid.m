function create_eigenvalues_grid()

   grid_eig_file = 'EigenvaluesInterpolation/eigenvalues_grid_low.csv';

   y0=0.025; ymax=0.3;
   for x = 1.05:1.05
	  
   for y = y0:0.005:ymax
     h = sqrt(y/y0) *0.0015;
     eig_values = get_sharp_eig_value(x,y,h);
     fid = fopen(grid_eig_file, "a");
     fprintf(fid, '%.3f, %.3f', x,y);
     fprintf(fid, ', %.17g', eig_values(1:4)');
     fprintf(fid, '\n');
     fclose(fid);
   end
   end

end
