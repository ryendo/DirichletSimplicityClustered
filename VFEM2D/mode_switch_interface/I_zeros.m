function A=I_zeros(m,n)
   global INTERVAL_MODE;
   if INTERVAL_MODE
      A = intval(zeros(m,n));
   else
      A = zeros(m,n);
   end
end



