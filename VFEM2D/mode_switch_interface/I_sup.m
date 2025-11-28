function output = I_sup(var)
   global INTERVAL_MODE;
   if INTERVAL_MODE
       if isintval(var)
          output = sup(var);
          return
       end
   end
  output = var;
end



