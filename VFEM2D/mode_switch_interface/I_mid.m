function output = I_mid(var)
   global INTERVAL_MODE;
   if INTERVAL_MODE
       if isintval(var)
          output = mid(var);
          return
       end
   end
  output = var;
end



