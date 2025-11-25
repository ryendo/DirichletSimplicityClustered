function A=I_sparse(m,n)
   global INTERVAL_MODE;
   if INTERVAL_MODE
      A = intval( sparse(m,n));
   else
      A = sparse(m,n);
   end
end



