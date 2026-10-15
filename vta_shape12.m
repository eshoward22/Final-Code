% ----- Spread, sphericity and Dice -----


function R = vta_shape12(Em, G, Eth, Em_ref)
%   Em_ref (optional) is the reference field (e.g. isotropic) for the Dice coefficient.

vta = (Em >= Eth);                 % NaN >= Eth is false, so the lead is never in the VTA
R.volume = vta_volume11(Em, G.h, Eth);


% ----- Spread -----
% The size of the VTA along x, y and z
[X, Y, Z] = ndgrid(G.x, G.y, G.z);
if any(vta(:))
    R.spread = [max(X(vta)) - min(X(vta)), max(Y(vta)) - min(Y(vta)), max(Z(vta)) - min(Z(vta))];
else
    R.spread = [0 0 0];
end

% ----- Sphericity -----
% Psi = pi^(1/3) (6V)^(2/3) / A_s
Ef = Em;
[Xr, Yr] = ndgrid(G.x, G.y);
rr = sqrt(Xr.^2 + Yr.^2);
a = max(rr(G.lead(:,:,end)));                  % lead radius on this grid
ring = (rr > a) & (rr <= a + 2.5*G.h);
for k = 1:length(G.z)
    slice = Ef(:,:,k);
    slice(G.lead(:,:,k)) = max(slice(ring));
    Ef(:,:,k) = slice;
end
Vencl = nnz(Ef >= Eth) * G.h^3;                % volume including the lead
Ef(isnan(Ef)) = 0;
% surface area of the |E| = Eth surface: isosurface() makes triangles, add their areas
[Xm, Ym, Zm] = meshgrid(G.x, G.y, G.z);        % isosurface needs meshgrid order, hence permute
fv = isosurface(Xm, Ym, Zm, permute(Ef, [2 1 3]), Eth);
if isempty(fv.faces)
    R.area = NaN;  R.sphericity = NaN;
else
    p1 = fv.vertices(fv.faces(:,1),:);
    p2 = fv.vertices(fv.faces(:,2),:);
    p3 = fv.vertices(fv.faces(:,3),:);
    tri = cross(p2 - p1, p3 - p1, 2);          % |cross product| / 2 = triangle area
    R.area = 0.5*sum(sqrt(sum(tri.^2, 2)));
    R.sphericity = pi^(1/3) * (6*Vencl)^(2/3) / R.area;
end

% ----- Dice coefficient -----
if nargin > 3
    ref = (Em_ref >= Eth);
    R.dice = 2*nnz(vta & ref) / (nnz(vta) + nnz(ref));
end
end
