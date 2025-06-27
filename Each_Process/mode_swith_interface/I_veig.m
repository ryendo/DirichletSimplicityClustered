function output=I_veig(A,B,indices)
   global INTERVAL_MODE;
   if INTERVAL_MODE
      [output,~] = veig(A,B,indices);
   else
      [~,d]=eigs(sparse(A),sparse(B),max(indices),'sm');
      output = diag(d);
   end
end



