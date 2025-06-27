function output=I_mid(I)
   global INTERVAL_MODE;
   if INTERVAL_MODE
      output = mid(I);
   else
      output = I;
   end
end



