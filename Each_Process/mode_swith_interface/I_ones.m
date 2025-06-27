function A=I_ones(m,n)
   global INTERVAL_MODE;
   if INTERVAL_MODE
      A = intval( ones(m,n));
   else
      A = ones(m,n);
   end
end



