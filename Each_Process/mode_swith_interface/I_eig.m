function output=I_eig(A,B,num_eigs)
   global INTERVAL_MODE;
   if INTERVAL_MODE

      %Use veigs for large matrix.
      [output,~] = veig(A,B,1:num_eigs);
   else
      [~,d]=eigs(sparse(A),sparse(B),num_eigs,'sm');
      output = diag(d);
   end
end



