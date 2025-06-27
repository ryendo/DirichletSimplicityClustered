function output = I_inf(var)
   global INTERVAL_MODE;
   if INTERVAL_MODE
      output = inf(var);
   else
      output = var;
   end
end



