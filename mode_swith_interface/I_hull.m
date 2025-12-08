function output=I_hull(I1,I2)
   global INTERVAL_MODE;
   if INTERVAL_MODE
      output = hull(I1,I2);
   else
      output = (I1+I2)/2;
      % output = I2;
      % output = (I1+I2)/2;
   end
end



