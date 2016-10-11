function [xbar,y,z,res] = admm_reformulated(c,A,b)
%% Solve a standard form LP using ADMM and reformulation techniques 
%% proposed by He & Yuan

maxiter = 10000; 
sigma = 10;
stoptol = 1e-5;
tau = 1.618; 
AT = A';
normb = norm(b);
normc = norm(c); 

[m,n] = size(A);
X = zeros(n,m);   % X = [x(1), ..., x(m)] 
Lam = zeros(n,m); % the lagrange multiplier
xbar = zeros(n,1); 
ee= ones(m,1); 
aa = sum(AT.*AT)';
breakyes = 0; 
for iter = 1:maxiter
    % Step 1
    invsigma = 1/sigma;
    Lamold = Lam;    
    H = (xbar-invsigma*c)*ee' - invsigma*Lam;    
    qtmp = (b - sum(AT.*H)')./aa;
    X = H + AT*spdiags(qtmp,0,m,m);
    % Step 2
    xbarinput = (X*ee+invsigma*Lam*ee)/m;
    xbar = max(xbarinput,0);
    z = sigma*(xbar-xbarinput); 
    y = (sigma/m)*qtmp;
    % Step 3
    Lam = Lamold + tau*sigma*(X -xbar*ee');
    % Check the primal-dual residual 
    primfeas = norm(b-A*xbar)/(1+normb);
    dualfeas = norm(c-AT*y-z)/(1+normc); 
    res(iter) = max(primfeas,dualfeas);
    if (res(iter) < stoptol)
       breakyes = 1; 
    end
    if (rem(iter,100)==1)
       primobj = c'*xbar;
       dualobj = b'*y;
       fprintf('\n %4.0f %3.2e %3.2e| %- 7.6e %- 7.6e| %3.2e',...
       iter,primfeas,dualfeas,primobj,dualobj,sigma);  
    end
    if (breakyes); break; end
end
%%*******************************************************************


