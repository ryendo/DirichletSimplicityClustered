function output = I_inf(var)
   global INTERVAL_MODE;
   if INTERVAL_MODE
       if isintval(var)
          output = inf(var);
          return
       end
   end
   output = var;
end



