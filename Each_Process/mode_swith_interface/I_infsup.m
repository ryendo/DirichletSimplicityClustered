function output=I_infsup(inf,sup)
   global INTERVAL_MODE;
   if INTERVAL_MODE
      output = infsup(inf,sup);
   else
      output = (inf+sup)/2;
   end
end



