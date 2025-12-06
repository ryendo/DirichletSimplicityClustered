function output=I_eigs(A,B,num_eigs,sigma)
   global INTERVAL_MODE;
   if isa(num_eigs, 'intval'), num_eigs = mid(num_eigs); end

   if INTERVAL_MODE            
      [output,~] = veigs(A,B,num_eigs,sigma);
   else
      [~,d]=eigs(sparse(A),sparse(B),num_eigs,sigma);
      output = diag(d);
   end
end



