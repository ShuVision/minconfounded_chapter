% test on real local covariance matrices from a subset of USPS dataset
% use shift-invert heuristic for trace ratio optimization
% trace_ratio_opt is slower than generalized eigensolver but results
% in significant larger ratio

clear;

% real matrix computed from a subset of usps dataset
fprintf('Load USPS\n');
load usps;

k = 5;
fprintf('\nk = %d\n\n',k);

tol = 1.e-6;

%%----------Maximize trace ratio v.s. max generalized eigs
fprintf('Maximizing tr(V^TAV)/tr(V^TBV) (with shift-invert heuristic)... ');
tic; [V, rho] = trace_ratio_opt(A,B,k,'max','tol',tol,'use_shift_invert',1,'verbose',1); t = toc;
fprintf(' tr(V^TAV)/tr(V^TBV) = %f. Time: %f\n', rho, t);

fprintf('Solving for largest generalized eigenvalues of (A,B)... ');
tic; [V,D] = eigs(A,B,k,'LM'); t = toc;
[V,R] = qr(V,0);
fprintf(' tr(V^TAV)/tr(V^TBV) = %f. Time: %f\n', trace(V'*A*V)/trace(V'*B*V), t);
