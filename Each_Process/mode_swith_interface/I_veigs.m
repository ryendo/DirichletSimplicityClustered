function output=I_veigs(A,B,SIGMA,indices)
   global INTERVAL_MODE;
   if INTERVAL_MODE
      [output,~] = veigs(A,B,SIGMA,indices);
   else
      [~,d]=eigs(sparse(A),sparse(B),max(indices),'sm');
      output = diag(d)';
   end
end



