%%**********************************************************************
%% Standard form LP: min c'x s.t. Ax=b
%% Dual: min (-b'y) s.t. A'y + z = c, z >= 0
%% We assume rank(A) = m
%%*********************************************************************
    function [x,y,z,res] ...
             = semi_proximal_alm_lp_solver(c,A,b)
    maxiter = 5000;
    stoptol = 1e-5;
    [m,n] = size(A);
    x = zeros(n,1); y = zeros(m,1); z = zeros(n,1);
    tau = 1.9;
    AAT = A*A';
    if (nnz(AAT)/prod(size(AAT)) > 0.4)
       L.R = full(chol(AAT)); L.matfct = 'chol';    
    else
       [L.R,indef,L.perm] = chol(AAT,'vector'); 
       L.Rt  = L.R'; L.matfct = 'spchol';     
    end
    normb = norm(b);
    normc = norm(c);
    breakyes = 0; 
    Rp = b-mexMatvec(A,x,0);
    Ac = mexMatvec(A,c,0); 
    Az = mexMatvec(A,z,0); 
    relgap = inf;
    sigma = 1; 
    for iter = 1:maxiter
        % step (1a)
        xold = x; 
        invsigma = 1/sigma;     
        rhs = invsigma*Rp + Ac - Az;
        ybar = linsysolve(L,rhs); 
        % step (1b)
        ATybar = mexMatvec(A,ybar,1);        
        z = max(c-ATybar-x/sigma,0);  
        Az = mexMatvec(A,z,0);         
        % step (1c)
        recompute = 1;
        if (recompute)
           rhs = invsigma*Rp + Ac - Az; 
           y = linsysolve(L,rhs); 
           ATy = mexMatvec(A,y,1);  
        else
           ATy = ATybar; y = ybar;
        end
        % step(2)
        Rd = ATy + z - c; 
        x = xold + tau*sigma*Rd;    
        
        % Calculate the current objective value
        Ax = mexMatvec(A,x,0); 
        Rp = b-Ax;
        primfeas = mexFnorm(Rp)/(1+normb);
        dualfeas = mexFnorm(Rd)/(1+normc);
        relcomp = abs(x'*z)/(1+norm(x)+norm(z));
        res(iter) = max(primfeas,dualfeas);
        if (res(iter) < stoptol) && (relcomp < stoptol) 
           breakyes = 1; 
        end;
        if (rem(iter,50)==1) || (breakyes)
           primobj = c'*x; dualobj = b'*y;
           relgap = abs(primobj-dualobj)/(1+abs(primobj));
           fprintf('\n %4.0f %3.2e %3.2e %3.2e| %- 7.6e %- 7.6e %3.2e|%3.2e',...
           iter,primfeas,dualfeas,relcomp,primobj,dualobj,relgap,sigma);
        end
        if (breakyes); break; end
        if (rem(iter,sigma_update_fun(iter))==0)
           primfeas2 = max(primfeas,0.1*stoptol);
           dualfeas2 = max(dualfeas,0.1*stoptol);
           ratio = primfeas2/dualfeas2;
           const = 1.1;
           if (max(ratio,1/ratio) > 500)
              const = const^3;
           elseif (max(ratio,1/ratio) > 50)
              const = const^2;
           end
           if (ratio > 5)
              sigma = max(1e-3,sigma/const);
           elseif (ratio < 1/5)
              sigma = min(1e3,sigma*const);
           end
        end
    end
 %%********************************************************************

   function sigma_update_iter = sigma_update_fun(iter)           
           
   const = 0.5;
   if (iter <= 25)
      sigma_update_iter = 10*const;
   elseif (iter <= 50)
      sigma_update_iter = 20*const;
   elseif (iter <= 100)
      sigma_update_iter = 40*const;
   elseif (iter <= 500)
      sigma_update_iter = 60*const;
   elseif (iter <= 1000)
      sigma_update_iter = 80*const;
   else
      sigma_update_iter = 100; 
   end
 %%********************************************************************



