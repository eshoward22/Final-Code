% ----- E = -grad V -----

function [E, G] = field_on_grid10(M, V, h, box)
%   Step 1: E at the centre of every element (from the shape functions).
%   Step 2: interpolate it onto a uniform grid with spacing h inside
%           box = [xmin xmax ymin ymax zmin zmax].
%   E.Ex, E.Ey, E.Ez, E.Em (= |E|) are arrays on the grid.
%   Grid points inside the lead are NaN (G.lead marks them).

% ----- Step 1 -----
% E at the element centres
X = M.X(M.conn);  Y = M.Y(M.conn);  Z = M.Z(M.conn);
U = V(M.conn);                                     % potential at the 8 corners
[Gx, Gy, Gz] = shape_gradients7(X, Y, Z, [0.5 0.5 0.5]);   % centre of the reference cube
Fc.Ex = reshape(-sum(Gx.*U, 2), M.ne);             % E = -sum_a V_a grad N_a
Fc.Ey = reshape(-sum(Gy.*U, 2), M.ne);
Fc.Ez = reshape(-sum(Gz.*U, 2), M.ne);
Fc.Em = sqrt(Fc.Ex.^2 + Fc.Ey.^2 + Fc.Ez.^2);
% The centres form a regular grid in (r, phi, z), so we can interpolate in those coordinates.

% ----- Step 2 -----
% The uniform sampling grid
G.x = box(1):h:box(2);  G.y = box(3):h:box(4);  G.z = box(5):h:box(6);  G.h = h;
[Xg, Yg, Zg] = ndgrid(G.x, G.y, G.z);
rg   = sqrt(Xg.^2 + Yg.^2);
phig = mod(atan2(Yg, Xg), 2*pi);
G.lead = rg < M.P.a;                                % inside the lead: not tissue
rq = max(rg, M.rc(1));                              % points between r = a and the first centre use the first centre

% phi is periodic, so add one extra layer of centres on each side (copied from the other end)
pc = [M.pc(end) - 2*pi, M.pc, M.pc(1) + 2*pi];
names = {'Ex', 'Ey', 'Ez', 'Em'};
for c = 1:4
    F = Fc.(names{c});
    F = cat(2, F(:,end,:), F, F(:,1,:));
    v = interpn(M.rc, pc, M.zcen, F, rq, phig, Zg, 'linear');   % trilinear interpolation
    v(G.lead) = NaN;
    E.(names{c}) = v;
end
end
