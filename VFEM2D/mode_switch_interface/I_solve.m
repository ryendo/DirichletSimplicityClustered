function output = I_solve(A,b)
   global INTERVAL_MODE;
   if INTERVAL_MODE
       output = full(A)\b;
   else
       output = A\b;
   end
end



