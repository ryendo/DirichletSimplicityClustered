function A=I_eye(m,n)
   global INTERVAL_MODE;
   if isa(m, 'intval'), m = mid(m); end
   if isa(n, 'intval'), n = mid(n); end
   if INTERVAL_MODE
      A = intval( eye(m,n));
   else
      A = eye(m,n);
   end
end



