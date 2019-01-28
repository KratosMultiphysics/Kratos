function info = trytime(A);
% TRYTIME : time SUPERLU with 3 outputs against Matlab's LU
%
% info = trytime(A);
% normally info is ||P*A*x-L*U*x|| / ||A|| (for a random x);
% but info is at least 10^6 if the factors are not triangular,
% or the permutation isn't a permutation.
% Copyright (c) 1995 by Xerox Corporation.  All rights reserved.
% HELP COPYRIGHT for complete copyright and licensing notice.


n = max(size(A));
info = 0;

format compact
disp('SUPERLU with 3 outputs:');

tic;[L,U,P]=lu(A);t1=toc;
fprintf('LU time = %d seconds.\n',t1);
tic;[l,u,pr]=superlu(A);t2=toc;
fprintf('SUPERLU time = %d seconds.\n',t2);
ratio = t1/t2;
fprintf('Ratio = %d\n',ratio);

if any(any(triu(l,1)))
    disp('L is *NOT* lower triangular.');
    info = info + 10^6;
else
    disp('L is lower triangular.');
end;
if nnz(l) == nnz(l+l)
    disp('L has no explicit zeros.');
else
    disp('L contains explicit zeros.');
    info = info+10^6;
end;
if any(any(tril(u,-1)))
    disp('U is *NOT* upper triangular.');
    info = info + 10^6;
else
    disp('U is upper triangular.');
end;
if nnz(u) == nnz(u+u)
    disp('U has no explicit zeros.');
else
    disp('U contains explicit zeros.');
    info = info+10^6;
end;
if pr == [1:n]
    disp('PROW is the identity permutation.');
elseif isperm(pr)
    disp('PROW is a non-identity permutation.');
else 
    disp('PROW is *NOT* a permutation.');
    info = info + 10^6;
end;

x = rand(n,1);
rnorm = norm(A(pr,:)*x - l*(u*x),inf)/norm(A,inf);
fprintf(1,'||A(PROW,:)*x -  L*U*x||/||A|| = %d\n', rnorm);

info = info + rnorm;
disp(' ');
