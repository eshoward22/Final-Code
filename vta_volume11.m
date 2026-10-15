% ----- VTA Volume -----

function vol = vta_volume11(Em, h, Eth)
%   (method breakdown, Section 4, Step 4).
%   Plain counting (number of voxels * h^3) is only O(h) accurate.  Instead:
%     d = (|E| - Eth) / |grad |E||        signed distance to the VTA boundary
%     f = min(max(1/2 + d/h, 0), 1)       fraction of the voxel that is activated
%     volume = sum(f) * h^3
%   NaN values (inside the lead) are not counted.
E = Em;
lead = isnan(E);
E(lead) = Eth;                  % placeholder so gradient() does not produce NaNs
E(isinf(E)) = 1e3;
[gy, gx, gz] = gradient(E, h);  % gradient() gives the derivative along the 2nd index first
g = sqrt(gx.^2 + gy.^2 + gz.^2);
d = (E - Eth) ./ max(g, eps);
f = min(max(0.5 + d/h, 0), 1);
f(lead) = 0;
vol = sum(f(:)) * h^3;
end
