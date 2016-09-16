density = @(M) nnz(A)/numel(A);

N=447610; 
A = sprandn(N,N,1/N)+sprandn(N,N,1/N)+sprand(N,N,1/N);
b = randn(N,1);
bnew = b+[zeros(N-100,1); randn(100,1)];

tic; 
x1 = A\b; x1new = A\bnew; 
toc;

% Using LUPQR
tic;
[L,U,P,Q,R] = lu(A); x2 = Q*(U\(L\(P*(R\b)))); x2new = Q*(U\(L\(P*(R\bnew))));
toc;
% Using LUPQ
tic;
[L,U,P,Q] = lu(A); x3 = Q*(U\(L\(P*b))); x3new = Q*(U\(L\(P*bnew)));

toc;
disp(norm(x1-x2)/norm(x1));
disp(norm(x1-x3)/norm(x1));
disp(norm(x1new-x2new)/norm(x1new));
disp(norm(x1new-x2new)/norm(x1new));
