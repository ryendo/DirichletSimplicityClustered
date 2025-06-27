function [neg_n,pos_n,zero_n,F]=GetInertia(D, AbortOnErr)
% Compute the inertia of diagonal matrix D.
%
% D is real symmetric 1*1 or 2*2 diagonal matrix.
%
%  In case the inertia can not be dertermined, if AbortOnErr==0, 
%  then F=1 will be returned, otherwise, the program exit with error.
%
% Xuefeng Liu, xfliu.math@gmail.com
% 2011/03/11 First version 
% 2011/11/22 Add parameter "AbortOnErr".

neg_n=0;
zero_n=0;
pos_n=0;

if( ~exist('AbortOnErr'))
    AbortOnErr = 0;
end

i=1;

while( i<= (size(D,1)-1) )

    if( D(i,i+1) == 0)
       if( D(i,i)>0) pos_n=pos_n+1;
       end
       if( D(i,i)<0) neg_n=neg_n+1;
       end
       i=i+1;
    else
        my_det=intval(D(i,i))*intval(D(i+1,i+1)) - intval( D(i,i+1))*intval( D(i+1,i));
        if( sup(my_det) <0 ) % two eigenvalue with different symbols.
            pos_n=pos_n+1;
            neg_n=neg_n+1;
            i=i+2; continue;
        end
        if( inf(my_det) > 0 ) % two eigenvalue with same symbols and non-zero.
            cond = intval( D(i,i) ) + intval( D(i+1,i+1) );
            if( inf(cond) >0 ) 
                pos_n=pos_n+2;
                i=i+2; continue;
            end
            if( sup(cond) <0 )
                neg_n=neg_n+2;
                i=i+2; continue;
            end
        end
        if(AbortOnErr==1)
            error('Can not determine the inertia of D ')
        else
            F=1
            return 
        end
    end
end

if( i == size(D,1) )
   if( D(i,i)>0) pos_n=pos_n+1;
   end
   if( D(i,i)<0) neg_n=neg_n+1;
   end
end
        
zero_n=size(D,1) - pos_n - neg_n;
F=0;
end
