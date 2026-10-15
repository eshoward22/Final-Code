% ----- 2D Finite Differences -----

function [V, Ex, Ey, Em, G, A] = fd_solve3(h, R, a, sx, sy, bc)
% FD_SOLVE  2D finite differences for   sx*V_xx + sy*V_yy = 0
%   bc(x,y) gives the potential at a boundary point, e.g. V0 on the
%   electrode and 0 on the outer circle.
%
%   Outputs:  V          potential at every grid point
%             Ex, Ey, Em electric field E = -grad V and its size |E|
%             G          the grid (x, X, Y, r, h, R, a, dom)
%             A          the finite-difference matrix (for the condition number)

% ----- Step 1 -----
% The grid with unknown points
x = (-R:h:R)';
N = length(x);
[X, Y] = ndgrid(x, x);            % X(i,j) = x(i),  Y(i,j) = x(j)
r = sqrt(X.^2 + Y.^2);
dom = (r > a) & (r < R);          % unknown points: between the electrode and the outer circle
n = nnz(dom);
num = zeros(N, N);
num(dom) = 1:n;                   % num(i,j) = equation number of grid point (i,j)
[I, J] = find(dom);               % grid position (i,j) of each unknown

% ----- Step 2 -----
% One equation per unknown point: five-point stencil and 
% the four neighbours: 1 = east, 2 = west, 3 = north, 4 = south
di  = [1 -1 0 0];
dj  = [0 0 1 -1];
sig = [sx sx sy sy];              % conductivity in that direction
opp = [2 1 4 3];                  % opposite neighbour (east <-> west, north <-> south)

rows = zeros(5*n, 1);  cols = rows;  vals = rows;   % entries of the sparse matrix A
b = zeros(n, 1);
k = 0;                                              % number of entries so far

for p = 1:n
    i = I(p);  j = J(p);
    xp = x(i); yp = x(j);

    % (a) distance to each neighbour. Normally h, but if the neighbour is past
    %     a boundary we use the point where the grid line crosses the circle.
    d   = [h h h h];
    isB = [false false false false];   % true if that "neighbour" is a boundary point
    Vb  = [0 0 0 0];                   % its (known) potential
    for m = 1:4
        if ~dom(i + di(m), j + dj(m))
            if r(i + di(m), j + dj(m)) <= a
                rho = a;               % the neighbour is inside the electrode
            else
                rho = R;               % the neighbour is outside the outer circle
            end
           
            B = xp*di(m) + yp*dj(m);
            C = xp^2 + yp^2 - rho^2;
            t = [-B - sqrt(B^2 - C), -B + sqrt(B^2 - C)];
            d(m)   = min(t(t > 0));    % first crossing in that direction
            isB(m) = true;
            Vb(m)  = bc(xp + d(m)*di(m), yp + d(m)*dj(m));
        end
    end

    % (b) the unequal-spacing formula (5) of the method breakdown:
    %     V_xx ~ 2/(hE+hW) * [ (V_E - V_P)/hE - (V_P - V_W)/hW ]
    %     so neighbour m gets weight 2*sigma/(d_m*(d_m + d_opposite)).
    %     When every d = h this is the usual five-point stencil.
    diag_val = 0;
    for m = 1:4
        w = 2*sig(m) / (d(m) * (d(m) + d(opp(m))));
        diag_val = diag_val - w;
        if isB(m)
            b(p) = b(p) - w*Vb(m);            % known value -> move to right-hand side
        else
            k = k + 1;
            rows(k) = p;
            cols(k) = num(i + di(m), j + dj(m));
            vals(k) = w;
        end
    end
    k = k + 1;
    rows(k) = p;  cols(k) = p;  vals(k) = diag_val;
end

A = sparse(rows(1:k), cols(1:k), vals(1:k), n, n);

% ----- Step 3 -----
% Solve A v = b (sparse direct solver)
u = A \ b;

% ----- Step 4 -----
% Put the answer back on the full grid
V = zeros(N, N);
V(r <= a) = bc(X(r <= a), Y(r <= a));    % electrode
V(r >= R) = bc(X(r >= R), Y(r >= R));    % outside the outer circle
V(dom) = u;

% ----- Step 1 -----
% Electric field E = -grad V (central differences)
% gradient() uses central differences inside and one-sided ones at the edges.
% It returns the derivative along the columns (here y) first.
[dVdy, dVdx] = gradient(V, h);
Ex = -dVdx;
Ey = -dVdy;
Em = sqrt(Ex.^2 + Ey.^2);

G.x = x; G.X = X; G.Y = Y; G.r = r; G.h = h; G.R = R; G.a = a; G.dom = dom;
end
