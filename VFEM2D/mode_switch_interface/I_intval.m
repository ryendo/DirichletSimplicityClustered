function output = I_intval(var)
   global INTERVAL_MODE;
   if INTERVAL_MODE
      output = intval(var);
   else
      output = var;
   end
end



