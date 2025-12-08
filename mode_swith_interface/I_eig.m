function output=I_eig(A,B,num_eigs)
   global INTERVAL_MODE;
   if isa(num_eigs, 'intval'), num_eigs = mid(num_eigs); end

   if INTERVAL_MODE
      [output,~] = veig(A,B,1:num_eigs);
   else
      [~,d]=eigs(sparse(A),sparse(B),num_eigs,'sm');
      output = diag(d);
   end
end



