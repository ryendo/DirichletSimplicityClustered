function A=I_sparse(m,n)
   global INTERVAL_MODE;
   if isa(m, 'intval'), m = mid(m); end
   if isa(n, 'intval'), n = mid(n); end
   if INTERVAL_MODE
      A = intval( sparse(m,n));
   else
      A = sparse(m,n);
   end
end



