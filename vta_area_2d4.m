% ----- 2D VTA Area -----


function [A_ray, A_vox, rth] = vta_area_2d4(Em, G, Eth, M)
% VTA_AREA_2D  Area of the 2D activated region {|E| >= Eth}, two ways.
%   (1) A_vox: count the activated grid points and multiply by h^2.
%   (2) A_ray: along M rays from the centre, find the activation radius
%       r_th(phi) by bisection, then  A = 1/2 * integral r_th^2 dphi
%       with the trapezium rule (method breakdown, Section 4, Step 5).
%   rth is the activation radius along each ray.
%   Note: A_ray includes the electrode disc; A_vox does not.
if nargin < 4
    M = 256;
end
h = G.h;

% (1) voxel count
active = (Em >= Eth) & (G.r > G.a) & (G.r < G.R);
A_vox = nnz(active) * h^2;

% (2) rays.  All M rays are bisected at the same time (lo, hi are vectors).
phi = (0:M-1)' * 2*pi/M;
f = @(rr) interpn(G.x, G.x, Em, rr.*cos(phi), rr.*sin(phi)) - Eth;   % |E| - Eth along the rays
lo = (G.a + 3*h) * ones(M, 1);      % near the electrode: |E| > Eth  (f > 0)
hi = (G.R - 3*h) * ones(M, 1);      % near the outer edge: |E| < Eth  (f < 0)
flo = f(lo);
for it = 1:45                        % 45 halvings: error <= (hi-lo)/2^46
    mid = (lo + hi)/2;
    fmid = f(mid);
    left = sign(fmid) == sign(flo);  % root is in the right half
    lo(left)  = mid(left);
    flo(left) = fmid(left);
    hi(~left) = mid(~left);
end
rth = (lo + hi)/2;

A_ray = 0.5 * sum(rth.^2) * (2*pi/M);   % trapezium rule (periodic)
end
