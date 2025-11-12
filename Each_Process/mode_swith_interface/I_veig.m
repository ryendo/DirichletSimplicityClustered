function output=I_veig(A,B,indices)
   global INTERVAL_MODE;
   if INTERVAL_MODE
      % Worse usage: 
      % [output,~] = veig(A,B,indices);

      %Use veigs for large matrix.
      [output,~] = veigs(A,B,'sm',max(indices));
   else
      [~,d]=eigs(sparse(A),sparse(B),max(indices),'sm');
      output = diag(d);
   end
end



