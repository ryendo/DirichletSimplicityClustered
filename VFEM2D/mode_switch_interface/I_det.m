function value = I_det(A)
   global INTERVAL_MODE;
   if INTERVAL_MODE
       n = size(A,1);
       if n > 3
           error('interval determinant only supports matrix with size <=3.')
       end
       A = intval(A);
       if n == 3
           tmp_1 = A(1,1)*A(2,2)*A(3,3) + A(1,2)*A(2,3)*A(3,1) + A(2,1)*A(3,2)*A(1,3);
           tmp_2 = A(1,3)*A(2,2)*A(3,1) + A(1,2)*A(2,1)*A(3,3) + A(1,1)*A(2,3)*A(3,2);
           value = tmp_1 - tmp_2;
           return
       end
       if n == 2
           value = A(1,1)*A(2,2)-A(1,2)*A(2,1);
       end
       if n ==1 
           value = A(1,1);
       end
   else
      value = det(A);
   end
end
