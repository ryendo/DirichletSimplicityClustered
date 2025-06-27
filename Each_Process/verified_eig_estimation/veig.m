function [eig_bounds, ind_range] = veig( A, B, ind )
%%  VEIG: Compute eigenvalue bounds for problem Ax=lambda Bx.
%  Input: 
%     A, B: symmetric real or interval matrix with size n;  B: positive definite. 
%     ind: (optional) the indices of eigenvalues wanted; default: 1:n.
%  Output:
%     eig_bounds: eigenvalue bounds in interval type;
%     ind_range : indices corresponding to returned eigenvalues.
%
%  Notice:
%  1) The returned result may have more eigenvalues than the ones specified by input parameter "ind".
%  2) This function is not expected to solve large scaled problem even the problem is a sparse one.
%  Reason: The block ldl faction method has been used to verify the bound of eigenvalues.
%  Although the block LDL algorithm can solve sparse problem efficiently, 
%  the calling of standard EIG function to compute ALL approximate eigenvalues
%  may take long time.
%
%  Reference: 
%  1) Yamamoto N., "A simple method for error bounds of eigenvalues of 
%  symmetric matrices", Linear Algebra and its Application, 234(2001), 227-234
%  2) Liu Xuefeng, Ogita takeshi, "An improvement of bounding eigenvalue of 
%  generalized eigenvalue problem", in preparation.
%
%  Xuefeng Liu, xfliu.math@gmail.com
%
%  2011/03/11 First version 
%  2011/09/30 Second version. Document modified.
%  2011/11/22 Third version. Modified to give exact bound if possible. Especially the part for repeated eigenvalues.
%
%%


%Check the symmetricity of matrix 
if( ~isequal(A,A'))
  error('Matrix A is not symmetric')
end
if( ~isequal(B,B'))
  error('Matrix B is not symmetric')
end


if(~ exist('ind') )
    min_ind=1; max_ind=size(A,1);
else
    if( ~isnumeric(ind) )
      error('The 3rd input parameter should be numeric.')
    end
    min_ind=min(ind); max_ind=max(ind);
end


eig_list = eig( full(mid(A)), full(mid(B)) );

eig_bounds = intval([]);

min_eig_B=min(eig(mid(B)));

index=min_ind;

while( index <= max_ind )

    
    eig_test = eig_list(index);

    
%%   Find rough lower bound and upper bound for index-th eigenvalue.


    lambda_test = eig_test;
    find_lambda = 0; k = 0;
    while( find_lambda < 1),
      if(k>34)
            error('Too many loops for finding rough lower bound')
      end
      E = mid(A) - lambda_test*mid(B);
        
      [L,D,P] = ldl(E);
      
      
      [neg_num,pos_num,zero_num,F] = GetInertia(D);      
      if( neg_num >= index || F~=0 )
        lambda_test = eig_test -max(eps,abs(eig_test))* 10^(k-16);
        k=k+1;
      else
        find_lambda = 1;
      end
    end

    neg_num_low = neg_num;
    lambda_test_low = lambda_test;

    
    
    lambda_test = eig_test;
    find_lambda = 0; k = 0;
    while( find_lambda < 1),
      if(k>34)
            error('Too many loops for finding rough upper bound')
      end
      E = mid(A) - lambda_test*mid(B);
      [L,D,P] = ldl(E);
      [neg_num,pos_num,zero_num,F] = GetInertia(D);      
      if( neg_num+zero_num < index || F~=0 )
        lambda_test = eig_test + max(eps,abs(eig_test)) * 10^(k-16);    k=k+1;
      else
        find_lambda = 1;
      end
    end


    neg_zero_num_upper = neg_num+zero_num;
    lambda_test_upper = lambda_test;

   
   %%%\-\-\-\-\\-\-\-\-\\-\-\-\-\\-\-\-\-\\-\-\-\-\\-\-\-\-\\-\-\-\-\\-\-%%%


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%             Find lower bound for index-th eigenvalue.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    lambda_test = lambda_test_low;
    E = mid(A) - lambda_test*mid(B);
    [L,D,P] = ldl(mid(E));
    E = intval(A) - lambda_test*intval(B);
    
    
    %%%%%%%%%%%%%%  Error estimation for lower bound %%%%%%%%%%%%%%

    diff_m = P*intval(L)*intval(D)*intval(L')*P' - E;

    diff_m_inf=norm(sup(abs(diff_m)),'inf');
    
    if( diff_m_inf ==0)
        lambda_lower = lambda_test;
        
    else
        err_est = max(eps,norm(sup(abs(diff_m)),'inf')/min_eig_B);
        err_est = max(eps, diff_m_inf /min_eig_B);
        is_positive=0;
        k=0;
        while( k<16 && is_positive<1)
            err_est = err_est*2;
            tmpMatrix = intval(err_est)*intval(B) - diff_m;	is_positive = isspd( hull(tmpMatrix , tmpMatrix') );
            tmpMatrix = intval(err_est)*intval(B) + diff_m;	is_positive = is_positive*isspd( hull(tmpMatrix , tmpMatrix') );
            clear tmpMatrix;
            k=k+1 ;
        end
        if( is_positive<1 )
            error('Too many loops for finding lower bound')
        end

        last_err_est = err_est;
        while( err_est > eps && is_positive)
            last_err_est = err_est;
            err_est = err_est/2;
            tmpMatrix = intval(err_est)*intval(B) - diff_m;	is_positive = isspd( hull(tmpMatrix , tmpMatrix') );
            tmpMatrix = intval(err_est)*intval(B) + diff_m;	is_positive = is_positive*isspd( hull(tmpMatrix , tmpMatrix') );
            clear tmpMatrix;
        end

        lambda_lower = inf( intval(lambda_test) - intval(last_err_est));
    end

    

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%             Find upper bound for index-th eigenvalue.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    lambda_test = lambda_test_upper;
    E = mid(A) - lambda_test*mid(B);
    [L,D,P] = ldl(mid(E));
    
    E = intval(A) - lambda_test*intval(B);


    %%%%%%%%%%%%%% Error estimate for upper bound %%%%%%%%%%%%%%

    diff_m = P*intval(L)*intval(D)*intval(L')* P' - E;

    diff_m_inf=norm(sup(abs(diff_m)),'inf');
    if( diff_m_inf ==0)
        lambda_upper = lambda_test;
    else
        err_est = max(eps, diff_m_inf /min_eig_B);

        is_positive=0;
        k=0;
        while( k<16 && is_positive<1)
            err_est = err_est*2;
            tmpMatrix = intval(err_est)*intval(B) - diff_m;	is_positive = isspd( hull(tmpMatrix , tmpMatrix') );
            tmpMatrix = intval(err_est)*intval(B) + diff_m;	is_positive = is_positive*isspd( hull(tmpMatrix , tmpMatrix') );
            clear tmpMatrix;
            k=k+1 ;
        end
        if( is_positive<1 )
            error('Too many loops for finding upper bound')
        end

        while( err_est > eps && is_positive)
            last_err_est = err_est;
            err_est = err_est/2;
            tmpMatrix = intval(err_est)*intval(B) - diff_m;	is_positive = isspd( hull(tmpMatrix , tmpMatrix') );
            tmpMatrix = intval(err_est)*intval(B) + diff_m;	is_positive = is_positive*isspd( hull(tmpMatrix , tmpMatrix') );
            clear tmpMatrix;
        end

        lambda_upper = sup( intval(lambda_test) + intval(last_err_est) );
    end
   %----------------------------------------------------------------------

    eig_bound = infsup(lambda_lower, lambda_upper);


    index_range=(neg_num_low+1):1:neg_zero_num_upper;
    
    if( neg_num_low+1 < min_ind ) 
        min_ind = neg_num_low+1;
    end
    
    eig_bounds( index_range-min_ind+1,1 ) = eig_bound;

    
    index = neg_zero_num_upper+1;

%%

end

ind_range = min_ind: (index-1); 

end

