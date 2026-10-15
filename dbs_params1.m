
% ----- Start Model Parameters -----
function P = dbs_params1(hc)

if nargin < 1
    hc = 0.25;
end

% ----- Electrode (Medtronic 3389) -----
P.a      = 0.635;        % lead radius (mm)
P.hc_len = 1.5;          % height of each contact (mm)
P.zc     = [0 2 4 6];    % centres of contacts c0, c1, c2, c3 (mm)

% ----- Tissue -----
P.sigma_iso = 0.2;       % average tissue conductivity (S/m)
P.sigma_enc = 0.1;       % encapsulation layer conductivity (S/m)
P.t_enc     = 0.5;       % encapsulation thickness (mm)
P.Eth       = 0.2;       % activation threshold (V/mm)

% ----- 3D mesh -----
P.hc     = hc;           % spacing in the core region (mm)
P.phi_r  = 4;            % radius (mm) where the arc length between angles equals hc
P.Rcore  = 8.5;          % radius of the evenly spaced core (mm)
P.Rout   = 40;           % radius of the grounded outer cylinder (mm)
P.zcore  = [-5 12];      % z-range of the evenly spaced core (mm)
P.zpad   = 35;           % extra length added above and below the core (mm)
P.growth = 1.3;          % each stretched layer is 1.3 times thicker than the last
end
