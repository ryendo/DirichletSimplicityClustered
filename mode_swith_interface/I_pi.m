function output = I_pi
   global INTERVAL_MODE;
   if INTERVAL_MODE
      output = intval(pi);
   else
      output = pi;
   end
end



