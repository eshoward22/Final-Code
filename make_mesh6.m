
% ----- 3D Mesh around the lead -----

function M = make_mesh6(P)

%   The mesh is built in cylindrical coordinates (r, phi, z):
%     - r:   from the lead surface r = a outwards, fine near the lead,
%            even in the core, stretched towards the outer boundary
%     - phi: Nphi equally spaced angles
%     - z:   even in the core (contact edges are mesh nodes), stretched beyond
%   Node (i,j,k) sits at (r_i cos(phi_j), r_i sin(phi_j), z_k).
%   Because the nodes lie exactly on r = a and r = a + t_enc, the round lead
%   and the encapsulation layer have no "staircase" edges.

a = P.a;  hc = P.hc;

% ----- Step 1 -----
% Radial nodes
r = a;
dr = hc/6;                                 % start very fine at the lead
while r(end) < a + P.t_enc - 1e-9          % (i) through the encapsulation layer,
    r(end+1) = min(r(end) + dr, a + P.t_enc);  %     growing 15% each step
    dr = dr*1.15;
end
while r(end) < P.Rcore - 1e-9              % (ii) keep growing until the spacing is hc
    dr = min(dr*1.15, hc);
    r(end+1) = min(r(end) + dr, P.Rcore);
end
dr = hc;
while r(end) < P.Rout - 1e-9               % (iii) stretch towards the outer boundary
    dr = dr*P.growth;
    r(end+1) = min(r(end) + dr, P.Rout);
    if P.Rout - r(end) < 0.5*dr            % avoid a tiny last element
        r(end) = P.Rout;
    end
end
M.rv = unique(r);

% ----- Step 2 -----
% Angles (a multiple of 8, so the x and y axes are mesh lines)
Nphi = 8*ceil(2*pi*P.phi_r/hc/8);
M.phiv = (0:Nphi-1) * 2*pi/Nphi;

% ----- Step 3 -----
% Z nodes
% (i) the core: split it into pieces at every contact edge, and put
%     about hc spacing inside each piece, so the contacts are exactly 1.5 mm
edges = unique([P.zcore(1), P.zc - P.hc_len/2, P.zc + P.hc_len/2, P.zcore(2)]);
edges = edges(edges >= P.zcore(1) & edges <= P.zcore(2));
z = edges(1);
for s = 1:length(edges)-1
    nseg = max(1, round((edges(s+1) - edges(s))/hc));
    piece = linspace(edges(s), edges(s+1), nseg + 1);
    z = [z, piece(2:end)];
end
% (ii) stretched layers below and above the core
below = [];  zl = z(1);  dz = hc;
while zl > P.zcore(1) - P.zpad + 1e-9
    dz = dz*P.growth;  zl = zl - dz;  below(end+1) = zl;
end
above = [];  zr = z(end);  dz = hc;
while zr < P.zcore(2) + P.zpad - 1e-9
    dz = dz*P.growth;  zr = zr + dz;  above(end+1) = zr;
end
M.zv = [fliplr(below), z, above];

% ----- Step 4 -----
% Node coordinates
nr = length(M.rv);  np = Nphi;  nz = length(M.zv);
[Rn, PHn, Zn] = ndgrid(M.rv, M.phiv, M.zv);
M.X = Rn(:).*cos(PHn(:));
M.Y = Rn(:).*sin(PHn(:));
M.Z = Zn(:);
% with this ordering, node (i,j,k) has number  i + (j-1)*nr + (k-1)*nr*np
node = @(i, j, k) i + (j-1)*nr + (k-1)*nr*np;

% ----- Step 5 -----
% Elements - each brick lists its 8 corner nodes
[I, J, K] = ndgrid(1:nr-1, 1:np, 1:nz-1);
I = I(:);  J = J(:);  K = K(:);
Jn = mod(J, np) + 1;                       % next angle (wraps from the last angle back to the first)
M.conn = [node(I,J,K),   node(I+1,J,K),   node(I,Jn,K),   node(I+1,Jn,K), ...
          node(I,J,K+1), node(I+1,J,K+1), node(I,Jn,K+1), node(I+1,Jn,K+1)];
% corner order: local node a = 1 + (step in r) + 2*(step in phi) + 4*(step in z)

% element centres (in r, phi, z) - used when sampling the field
M.rc = (M.rv(1:end-1) + M.rv(2:end))/2;
M.pc = M.phiv + pi/Nphi;
M.zcen = (M.zv(1:end-1) + M.zv(2:end))/2;
M.ne = [nr-1, np, nz-1];                   % number of elements in r, phi, z

% which elements are in the encapsulation layer
rc = M.rc(:);
M.enc = rc(I) < a + P.t_enc - 1e-9;

% ----- Step 6 -----
% Nodes where the potential will be fixed
Rall = Rn(:);  Zall = Zn(:);
M.boundary = find(Rall >= P.Rout - 1e-9 | Zall <= M.zv(1) + 1e-9 | Zall >= M.zv(end) - 1e-9);  % grounded
onLead = Rall <= a + 1e-9;
M.leadsurf = find(onLead);                 % the whole lead surface
M.contacts = cell(1, length(P.zc));
for c = 1:length(P.zc)                     % the four contacts
    M.contacts{c} = find(onLead & abs(Zall - P.zc(c)) <= P.hc_len/2 + 1e-9);
end
M.P = P;
end
