% ----- Boundary Conditions -----

function [V, info] = solve_fem9(K, fixed, values, tol)
%   fixed  - the nodes with a known potential (active contact, outer boundary)
%   values - the potentials at those nodes
%   The rest of the system is symmetric positive definite, so it is solved
%   with preconditioned conjugate gradients (PCG) and an incomplete
%   Cholesky preconditioner.
if nargin < 4
    tol = 1e-10;
end
N = size(K, 1);
V = zeros(N, 1);
V(fixed) = values;
free = true(N, 1);
free(fixed) = false;

% Split K V = 0 into known and unknown parts:  K_ff V_f = -K_fd V_d
Kff = K(free, free);
rhs = -K(free, fixed) * values(:);

t0 = tic;
% diagonal scaling: A = D K_ff D with D = diag(1/sqrt(K_ii)), so A has 1's on the diagonal
D = spdiags(1./sqrt(full(diag(Kff))), 0, nnz(free), nnz(free));
A = D*Kff*D;
b = D*rhs;
L = ichol(A, struct('type', 'ict', 'droptol', 1e-3));   % incomplete Cholesky: A ~ L L'
[y, flag, relres, iter] = pcg(A, b, tol, 3000, L, L');
V(free) = D*y;

info.flag = flag; info.relres = relres; info.iter = iter;
info.time = toc(t0); info.nfree = nnz(free);
end
