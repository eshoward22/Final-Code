% ----- Conductivity Tensor -----

function S = sigma_tensor2(sigma_iso, kappa, theta_deg)
% The 3x3 conductivity tensor 
%   kappa     = sigma_par / sigma_perp   (anisotropy ratio)
%   theta_deg = fibre angle in degrees, rotating about the x-axis
%               (theta = 0: fibres along the lead, theta = 90: fibres along -y)
%   kappa = 1 gives the isotropic tensor sigma_iso * I.

% ----- Step 1 -----
% The two conductivities, using the volume constraint: sigma_par * sigma_perp^2 = sigma_iso^3
s_perp = sigma_iso / kappa^(1/3);
s_par  = kappa * s_perp;

% ----- Step 2 -----
% The tensor when the axes line up with the fibres (fibres along z)
S_aligned = diag([s_perp, s_perp, s_par]);


% ----- Step 3 -----
% Rotate it to the fibre angle theta
th = theta_deg * pi/180;
Q = [1    0         0;
     0  cos(th)  -sin(th);
     0  sin(th)   cos(th)];
S = Q * S_aligned * Q';

S = (S + S')/2;   
end
