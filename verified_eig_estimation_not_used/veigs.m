function [lambda,ind_range] = veigs( A, B, SIGMA, NumOfEigs )
%% VIEGS(A,B,SIGMA,NumOfEigs) finds verified bounds for eigenvalues of Ax=lambda Bx.
%
% ------------- Input Parameters: A, B, SIGMA  ----------------------------
%
% A and B  must be symmetric real or interval matrices and B be positive 
% definite. If B is not positive definite, the program will fail to continue;
% A and B can be in sparse matrix format. The sheme for calculating approximate
% eigenvalue is dependent on format of A and B. 
%
% 
% SIGMA is optional. Default value of SIGMA is 'lm'.
% NumOfEigs, optional one, speifies the least number of eigenvalues to be 
% returned. Default value is 1. 
% VEIGS(A,B,SIGMA,NumOfEigs) returns bound for eigenvalue(s) as below:
% SIGMA = 'sa':  smallest eigenvalue; 
% SIGMA = 'la':  largest eigenvalue; 
% SIGMA = 'sm':  eigenvalue of smallest magnitude;
% SIGMA = 'lm':  eigenvalue of largest magnitude;
% If SIGMA is a real scalar, VEIGS(A,B,SIGMA) find the bound for eigenvalue(s) near to SIGMA. 
% ----------------------- Output ------------------------------------------
%
% [lambda] = VEIGS(A,B) return the bound(s) of eigenvalues in INTVAL
% format. The returned lambda can be one eigenvalue or several eigenvalue
% in a cluster.
%
% [lambda, ind] = VEIGS(A,B) returns both eigenvalue bounds and correspoing 
% index(indices).
% 
% If either A or B is in sparse matrix format, VEIGS first call EIGS to 
% compute the approximate eigenvalues. If EIGS failed or A and B are fullmatrix, 
% the code turns to EIG, which may take longer time for large matrix.
%
%
% Exmaple 1:
%     n=8; A=eye(n); B=infsup(hilb(n)-1E-13,hilb(n)+1E-13);
%     [bounds, ind]=veigs(A,B,0.1)
%     [bounds, ind]=veigs(A,B,'sa')
%     [bounds, ind]=veigs(A,B,'la')
% Example 2:
%     n=10; A=diag(ones(n,1))*6; A(1,1)=5; A(n,n)=5;
%     A=A + diag( ones(n-1,1),1)*(-4) + diag( ones(n-1,1),-1)*(-4);
%     A=A + diag( ones(n-2,1),2) + diag( ones(n-2,1),-2);
%     B=hilb(n)*232792560;
%     [bound, ind]=veigs(A,B,'la')
%
% Xuefeng Liu, xfliu.math@gmail.com
% First Version. March 4th, 2011
% 
% Reference: H. Behnke, Clausthal, The calculation of Guaranteed bounds for
% eigenvalues using complementary variational principles, Computing 47, 
% 11-27 (1991).
% (Multiple precision computation version sugguested in above paper has not been implemented.)
%
% Last updated: Oct. 5 2011 
% Last updated: Nov. 22 2011 : Modified to deal with repeated eigenvalues.
%
%


%% --------- Configuration of parameters -----------------------------------------

%Selection of method to prove the existence of eigenvalues.
% 1: method of Rump (command:verifyeig)
% 2: method of Behnke, where Lehmann's theorem is the keypoint.
MethodSelect = 1;

% The parameter lambda_ratio and do_shift are for Behnke's method.

% The value of lambda_ratio is used to determine rho and sigma.
% rho = eig_list(r-1)*(1-lambda_ratio) + eig_list(r)*lambda_ratio;
% sigma = eig_list(s)*lambda_ratio + eig_list(s+1)*(1-lambda_ratio);
% where eig_list is a local list of eigenvalues.
lambda_ratio = 1E-6;

% The minimum relative distance among eigenvalue in a cluster. 
% lambda_r - lambda_{r+1} >  min_dist ( lambda_r + lambda_{r+1})
min_dist = 1E-8;

% Wheather eigenvalue shift is used or not.
do_shift = 1; 


%%% ----------- Step 1: parsing parameter -------------- %%%%
lambda_shift = 0.0;
lambda = nan; ind_range = nan;

eig_list = []; %Local list of eigenvalue in increasing order.
SIGMA_max = 0; SIGMA_min = 0; 

if( ~exist('SIGMA') )
   SIGMA='lm';
end
% The number of approximate eigenvalues required in calling "eigs".
if( ~exist('NumOfEigs'))
   NumOfEigs=1;
   EigNum = min( 10, size(A,1) );  
else
   EigNum = min(size(A,1), NumOfEigs+2);
end


if( isstr(SIGMA) )
    if( strcmp( lower(SIGMA), 'sm')  || strcmp( lower(SIGMA), 'lm') || strcmp( lower(SIGMA), 'la')  || strcmp( lower(SIGMA), 'sa') )
        [V,D,F] = eigs( mid(A),mid(B),EigNum, SIGMA );   
        if( F ~= 0 ) %'Eigenvalue does not convergent. We turn to use eig() '
            [V,D] = eig( full(mid(A)),full( mid(B)));
        end
        [eig_list,p]=sort(diag(D));
    else
        error(['Parameter SIGMA not available: ', SIGMA]) ;
    end

    if( strcmp( lower(SIGMA), 'sa') )
        SIGMA_min = 1;
        [v,ind] =  min(eig_list);
    elseif( strcmp( lower(SIGMA), 'la') )
        SIGMA_max = 1;
        [v,ind] =  max(eig_list);
    elseif( strcmp( lower(SIGMA), 'sm') )
        [v,ind] =  min( abs(eig_list) );
    else % 'lm'
        [v,ind] =  max( abs(eig_list) );
    end
else
    [V,D,F] = eigs( mid(A),mid(B),EigNum, SIGMA );   
    if( F ~= 0 ) %'Eigenvalue does not convergent.'
        [V,D] = eig( full(mid(A)),full( mid(B)));
    end
    [eig_list,p]=sort(diag(D));
    [v,ind] =  min( abs(eig_list - ( SIGMA ) ));   
    
end

V=V(:, p);
m=size(V,2); % eig_list(m) has the max value of eig_list


global_sigma = inf; global_rho = -inf;
global_s = -1; global_r = size(A,1)+1;


%%% ----------- Step 2: determination of r and s -------------- %%%%
%% ------------------------------------------------------------------
% Find the index r,s.
% Refer to paper of Behnke for definition of r and s.
% lambda_{r-1} < lambda_{r} approx lambda_{s} < lambda_{s+1}
% The (ind)-th  eigenvalue is the approximate eigenvalue of interest.
% ------------------------------------------------------------------



%% ------------------------------------------------------------------
%                    Determination of s
%  ------------------------------------------------------------------
s=ind;

if(SIGMA_max ==1) % LM: largest eigenvalue.
    s=m;
    global_sigma = veigs_RoughUpper( A, B, eig_list(s),V(:,s), lambda_ratio);
    global_s=size(A,1);
end
if (SIGMA_min ==1 ) %SM: smallest eigenvalue;
    r=1;
    global_rho = veigs_RoughLower(A,B, eig_list(r),V(:,r), lambda_ratio);
    global_r=1;
end


if( global_sigma == inf )

    while s<m &&  (eig_list(s) > eig_list(s+1) - lambda_ratio*abs(eig_list(s+1)) || s-ind < NumOfEigs-1 )
        s=s+1;
    end

    if( s==m)
        eig_next = eig_list(s) + abs(eig_list(s))*0.01; 
    else
        eig_next = eig_list(s+1);
    end
    
    lambda = lambda_ratio*eig_list(s) + (1-lambda_ratio)*eig_next;

    tmpM = intval(A)-lambda*intval(B);

    if exist('ldl')
      [L,D,P]= ldl( mid(tmpM) );
    else
      [L,D,P]= my_ldl( mid(tmpM) );
    end
    [neg_num,pos_num,zero_num] = GetInertia(D,1);
    
    if( neg_num + zero_num == size(A,1) ) 
       SIGMA_max = 1;
       s=m; global_sigma = veigs_RoughUpper( A, B, eig_list(s),V(:,s), lambda_ratio); global_s=size(A,1);
    else

        diff_m = P*intval(L)*intval(D)*intval(L')*P' - tmpM;%The result matrix becomes not symmetric.

        err_est_base = max(eps, (1-lambda_ratio)*( eig_next - eig_list(s) )/2^4 ); 
        err_est = err_est_base;
        is_positive = 0;    scale = 1;
        while(  is_positive == 0 )
            err_est = err_est_base*(scale);
            if( err_est > (1-lambda_ratio)*( eig_next - eig_list(s) ) && s<m ),
                 break;
            end
            % Estimation for | lambda_i( E, B) - lambda_i( PLDL'P', B) | < err_est  
            %by confirming err_est*B - (E-PLDL'P')  and err_est*B - (PLDL'P'-E) are both positively definite;
            tmpMatrix = intval(err_est)*intval(B) - diff_m; is_positive = isspd( hull(tmpMatrix, tmpMatrix') );
            tmpMatrix = intval(err_est)*intval(B) + diff_m;	is_positive = is_positive*isspd( hull(tmpMatrix, tmpMatrix') );
            clear tmpMatrix ;
            scale = scale*2; 
        end
        if(~ is_positive )
            if( neg_num == 1 ) % We turn to a rough lower bound by a direct isspd test.
                %It arrived here because it failed to give the upper bound
                %for minimum eigenvalue by LDL. But we can get such bound
                %by Rayleigh bound, which will be done in "LDL Result check"
                %section.
                SIGMA_min = 1;
                r=1; global_rho = veigs_RoughLower(A,B, eig_list(r),V(:,1), lambda_ratio); global_r=1;
            else
            error( ['LDL fraction: failed to find lower bound of eigenvalue lambda_{',num2str(s+1),'}.'] );
            end
        else
            global_s = neg_num;
            global_sigma = inf( lambda - intval(err_est) ); 
            %
            % The (global_s + 1)th eigenvalue must be greater than global_sigma. 
            % That is to say, there exists at maximum (global_s) eigs in (-inf, global_sigma ]
            % Lehmann's theorem will declare the existence of eig in 
            % (global_sigma + 1/mu_i, global_sigma), thus the lower bounds for (global_s+1)th eigenvalue.
            %
        end
    end
end

%% -------------------------------------------------------------------------
% Determination of r
% --------------------------------------------------------------------------

if( global_rho == -inf )
    r=ind;
    while r>1 && ( eig_list(r) < eig_list(r-1) + lambda_ratio*abs(eig_list(r-1)) || ind -r < NumOfEigs-1 )
            r=r-1;
    end

    if r==1
        eig_pre = eig_list(r) - abs(eig_list(r))*0.01;
    else
        eig_pre = eig_list(r-1);
    end

    lambda = lambda_ratio*eig_list(r) + (1-lambda_ratio)*eig_pre;
    tmpM = intval(A)-lambda*intval(B);

    if exist('ldl')
      [L,D,P]= ldl( mid(tmpM) );
    else
      [L,Tmp,P]= chol( sparse(mid(tmpM)), 'lower' );
      if Tmp >0 
        error('Failed in chol computing.')
      end
      D = diag(diag(L).^2);
      L = L*diag( 1./diag(U) );
    end

    neg_num = GetInertia(D);

   
    if( neg_num == 0 )
       SIGMA_min = 1; 
       r=1; global_rho = veigs_RoughLower(A,B, eig_list(r),V(:,1), lambda_ratio); global_r=1;
    else
        diff_m = P*intval(L)*intval(D)*intval(L')*P' - tmpM;
        err_est_base = max(eps,(1-lambda_ratio) * ( eig_list(r) - eig_pre )/2^4); 
        err_est = err_est_base;
        is_positive = 0;    scale = 1;
        while(  is_positive == 0 )
            err_est = err_est_base*(scale);
            if( err_est > (1-lambda_ratio)*( eig_list(r) - eig_pre ) && r>1 )
                break;
            end
            % Estimation for | lambda_i( E, B) - lambda_i( PLDL'P', B) | < err_est  
            % by confirming err_est*B - (E-PLDL'P')  and err_est*B - (PLDL'P'-E) are both positively definite;
            tmpMatrix = intval(err_est)*intval(B) - diff_m;	is_positive = isspd( hull(tmpMatrix , tmpMatrix') );
            tmpMatrix = intval(err_est)*intval(B) + diff_m;	is_positive = is_positive*isspd( hull(tmpMatrix , tmpMatrix') );
            clear tmpMatrix;
            scale = scale*2;
        end

        if(~ is_positive) % We turn to a rough upper bound by a direct isspd test
            if( neg_num == size(A,1)-1 ) 
                %It arrived here because it failed to give the lower bound for the maximum eigenvalue
                SIGMA_max = 1;
                s=m; global_sigma = veigs_RoughUpper( A, B, eig_list(s),V(:,s), lambda_ratio); global_s=size(A,1);
            else
                error(['LDL fraction: failed to find upper bound of eigenvalue lambda_{',num2str(r-1),'}.']) ; 
            end
        else
            global_r = neg_num+1;
            global_rho = sup(lambda + intval(err_est)); %the selection of lambda and err_est_base make global_rho < eig(global_r).
        end

    end

end

%% -------------------------------------------------------------------------
%  LDL Result check
%  Check whether LDL failed or not.
%%-------------------------------------------------------------------------

if( global_s - global_r ~= s - r  || global_s < 0 || global_r > size(A,1) )
    if( SIGMA_max ==1 ) % Rayleigh bound for maximum eigenvalue
        lower_bound = V(:,m)'*intval(A)*V(:,m)/(V(:,m)'*intval(B)*V(:,m));
        lambda = infsup( lower_bound.inf, global_sigma);
        ind_range = size(A,1);
        return
    end
    if( SIGMA_min ==1 ) % Rayleigh bound for minimum eigenvalue
        
        upper_bound = V(:,1)'*intval(A)*V(:,1)/(V(:,1)'*intval(B)*V(:,1));
        lambda = infsup( global_rho, upper_bound.sup);
        ind_range = 1
        return
    end
end

if( global_s - global_r > s - r )
    error('Computation failed in finding upper bound and lower bound for cluster eigenvalues around given value. Please consider to increase value of EigNum in this code.')
end
if( global_s - global_r < s - r )  %This happens in very rare case when approximate eigenvalues are in poor precision
    error('Computation failed in finding rough  bounds  rho < lambda_s ~ lambda_r < sigma ')
end


%When we apply Lehmann's theorem, it is better to have sigma and rho departed from the object eigenvalue (cluster).
if( SIGMA_max ==1)
    global_sigma = max( global_sigma, eig_list(m) + abs(eig_list(m) ));
end
if( SIGMA_min ==1)
    global_rho = min( global_rho, eig_list(r) - abs( eig_list(r)) );
end

%%% global_sigma: lower bound of s+1 th eigenvalue .
%%% global_rho:  upper bound of r-1 th eigenvalue.


V=V(:,r:s);

%%% Lehmann-Behnke's method.
lambda = Lehmann_Behnke__(A,B, eig_list, V, global_rho, global_sigma, r,s, do_shift);

%%% Apply Rump's method ( verifyeig ) by checking the eigenvalue existence interval.
% lambda = Rump__(A,B,eig_list,V,global_rho, global_sigma,r,s);

ind_range = global_r:1:global_s;

end




%%  By checking isspd( A - eig_test * B) to find rough lower bound.
function [bound]=veigs_RoughLower(A,B,eig_test,v,lambda_ratio)

    k=1;
    delta = max(eps, abs(eig_test)/2^4);
    min_test = eig_test - max(1, abs(eig_test));
    lambda = eig_test - delta*2^(k-1);
    while( lambda > min_test || k==1)
        tmpMatrix =  intval(A) - lambda * intval(B) ;
        is_pd = isspd( hull(tmpMatrix,tmpMatrix') );
        if( is_pd) 
            break
        end
        k=k+1; lambda = eig_test - delta*2^(k-1);
    end
    if(is_pd)
        bound = lambda;
    else
        error('Failed to find lower bound for minimum eigenvalue'); 
    end


end

%% By checking isspd( eig_test * B - A ) to find rough upper bound.
function [bound] = veigs_RoughUpper(A,B,eig_test,v,lambda_ratio)
    k=1;
    delta = max(eps, abs(eig_test)/2^4);
    max_test = eig_test + max(1, abs(eig_test));
    lambda = eig_test+ delta*2^(k-1);
    while( lambda < max_test || k==1)
        tmpMatrix =  lambda * intval(B) - intval(A) ;
        is_pd = isspd( hull(tmpMatrix,tmpMatrix') );
        if( is_pd) 
            break
        end
        k=k+1; lambda = eig_test+ delta*2^(k-1);
    end
    if(is_pd)
        bound = lambda;
    else
        error('Failed to find upper bound for maximum eigenvalue'); 
    end

end


%% 
function lambda=Lehmann_Behnke__(A,B, eig_list, V, global_rho, global_sigma,r,s,do_shift)

%%

%  Lower bound for min eig(B);

[v_test, c_test,F] = eigs(mid(B),1,'sm');
if( F ~= 0 ) %'Eigenvalue does not convergent. We turn to use eig() '
    [v_test, c_list] = eig(mid(B));
    c_test = min(diag(c_list));
end
if( c_test <= 0)
    error('Second parameter matrix should be positive definite.');
end

clear v_test;

c_test=intval(c_test)*0.95;
while( ~isspd( intval(B) -  c_test*speye( size(B,1) )) && c_test > eps)
    c_test = c_test*0.5;
end

if (c_test  <= eps )
    error('Failed to find lower bound for eigvalues of inputed B'); 
end
c=c_test;

%%
%  ----------- Lehmann's bounds -------------- %



lambda_shift = 0;
if do_shift>0,
    lambda_shift = eig_list(r);
    A = intval(A) - lambda_shift*intval(B);
    global_rho = sup(global_rho-intval(lambda_shift));
    global_sigma = inf(global_sigma-intval(lambda_shift));

end


%------------  creating small submatrix needed by Lehmann's method -----

A0 = V'*intval(B)*V;  %interval
A1 = V'*intval(A)*V;  %interval
nV = mid(B)\(mid(A)*mid(V));
err_vec = intval(B)*intval(nV) - intval(A)*intval(V);
Err = err_vec'*err_vec/intval(c);

% Here, inverse B is needed. Since err_vec is small, a method that gives 
% rough estimate is also acceptable. For example: 
% Err = err_vec'*( intval(B) \ err_vec);


A2 = V'*A*intval(nV) - intval(nV)'*( B*intval(nV) - A*V ) + Err;


% Computation of lower bound

sigma = global_sigma;

SA = (- A2 + sigma*A1 );
SB = (- A1 + sigma*A0 );

if ~ isspd( hull(SB,SB') )
    error('Error in applying Lehmann theorem.');
end

if( s==r )
    low =   SA/SB ;
else
    low =   veig( hull(SA,SA'), hull(SB,SB'));
end


% Computation of upper bound
%upper bound: A-> -A; lambda-> -lambda


rho =  global_rho;

SA = A2-rho*A1;
SB = A1-rho*A0;

if ~ isspd( hull(SB,SB') )
    error('Error in applying Lehmann theorem.');
end

if( s == r  )
    upper = (SA/SB) ;    
else
    upper =   veig( hull(SA,SA'), hull(SB,SB'));; 
end

lambda = infsup( inf(low+intval(lambda_shift)), sup(upper+intval(lambda_shift)));


end


%% 
function lambda = Rump__(A,B,eig_list,V,global_rho, global_sigma,r,s)

if( r>s)
    approx_eig = 0.5*global_rho+0.5*global_sigma ;
else
    approx_eig = eig_list(r);
end
    
eig_bounds = verifyeig(A,approx_eig, V,B);
if( max( sup(eig_bounds)) < global_sigma &&   min( inf(eig_bounds)) > global_rho )
    lambda = eig_bounds;
else
    lambda = eig_bounds;
    display('The indices for eigevanlues can not be verified.');
end

end

