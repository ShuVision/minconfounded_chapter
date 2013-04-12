% test trace_ratio_opt on randomly generated matrices
% trace_ratio_opt is faster than generalize eigensolver
clear;

fprintf('Generating random symmetric positive definite matrices A and B\n');
n = 3000;
d = 20;
U1 = qr(rand(n,d),0); A = U1*diag([1:d].^3)*U1';
U2 = rand(n); B = U2*U2';
% sparsify the matrices, comment the below 2 lines to test on sparse matrices
%A = sparse((rand(n)>0.2).*A);
%B = sparse((rand(n)>0.2).*B);
% symmetrize the matrices
A = (A+A')/2;
B = (B+B')/2;
% regularize to make sure the matrices are positive
A = A + 0.1*diag(repmat(abs(trace(A)), n, 1)); 
B = B + 0.1*diag(repmat(abs(trace(B)), n, 1));

k = 10;
fprintf('\nk = %d\n\n',k);

tol = 1.e-6;
opts.tol = tol;

%%----------Maximize trace ratio v.s. max generalized eigs
fprintf('Maximizing tr(V^TAV)/tr(V^TBV)... ');
tic; [V, rho] = trace_ratio_opt(A,B,k,'max','tol',tol,'use_shift_invert',0); t = toc;
fprintf(' tr(V^TAV)/tr(V^TBV) = %f. Time: %f\n', rho, t);

fprintf('Solving for largest generalized eigenvalues of (A,B)... ');
tic; [V,D] = eigs(A,B,k,'LM',opts); t = toc;
[V,R] = qr(V,0);
fprintf(' tr(V^TAV)/tr(V^TBV) = %f. Time: %f\n', trace(V'*A*V)/trace(V'*B*V), t);

fprintf('\n');

%%----------Minimize trace ratio v.s. min generalized eigs
fprintf('Minimizing tr(V^TAV)/tr(V^TBV)... ');
tic; [V, rho] = trace_ratio_opt(A,B,k,'min','tol',tol,'use_shift_invert',0,'use_adaptive_tol',0); t = toc;
fprintf(' tr(V^TAV)/tr(V^TBV) = %f. Time: %f\n', rho, t);

fprintf('Solving for smallest generalized eigenvalues of (A,B)... ');
tic; [V,D] = eigs(A,B,k,'SM',opts); t = toc;
[V,R] = qr(V,0);
fprintf(' tr(V^TAV)/tr(V^TBV) = %f. Time: %f\n', trace(V'*A*V)/trace(V'*B*V), t);

