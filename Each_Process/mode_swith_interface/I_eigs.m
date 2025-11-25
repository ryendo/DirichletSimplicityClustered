function output=I_eigs(A,B,num_eigs,sigma)
   global INTERVAL_MODE;
   if INTERVAL_MODE

      %Use veigs for large matrix.
      [output,~] = veigs(A,B,num_eigs,sigma);
   else
      [~,d]=eigs(sparse(A),sparse(B),num_eigs,sigma);
      output = diag(d);
   end
end



