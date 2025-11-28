function A=I_eye(m,n)
   global INTERVAL_MODE;
   if INTERVAL_MODE
      A = intval( eye(m,n));
   else
      A = eye(m,n);
   end
end



