function output = I_intval(var)
   global INTERVAL_MODE;
   if INTERVAL_MODE
      output = intval(var);
   else
      if ischar(var)
        output = str2double(var);
      else
        output = var;
      end
   end
end



