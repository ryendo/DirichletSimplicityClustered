function output = I_inf(var)
   global INTERVAL_MODE;
   if INTERVAL_MODE
       if isintval(var)
          output = sup(var);
       else
          output = var;
       end
       
   else
      output = var;
   end
end



