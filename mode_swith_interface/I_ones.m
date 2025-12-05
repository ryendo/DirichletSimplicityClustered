function A=I_ones(m,n)
   global INTERVAL_MODE;
   if isa(m, 'intval'), m = mid(m); end
   if isa(n, 'intval'), n = mid(n); end
   
   if INTERVAL_MODE
      A = intval( ones(m,n));
   else
      A = ones(m,n);
   end
end



