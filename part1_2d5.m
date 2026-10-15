% ----- 2D model -----
%  Solves  sx*V_xx + sy*V_yy = 0  around a disc electrode (radius a, held at V0)
%  inside a grounded circle (radius R).  The fibres run along x.
%
%  A  Manufactured solution  -> is the solver really 2nd order?
%  B  Analytic isotropic VTA -> grid convergence of the VTA area + Richardson
%  C  Trapezium rule         -> how many rays are needed for the area?
%  D  Truncation vs round-off error in the central difference for E
%  E  Condition number of the linear system
%  F  Domain size            -> does the distance to the outer boundary matter?
%  G  Main 2D comparison     -> isotropic vs anisotropic (Figure 3)

clear; close all;
if ~exist('Finalresults', 'dir'), mkdir('Finalresults'); end
if exist('Finalresults/summary2d.txt', 'file'), delete('Finalresults/summary2d.txt'); end
diary('Finalresults/summary2d.txt');      % copies everything printed into this file

% ----- Parameters -----
V0  = 2;          % electrode voltage (V)
a   = 0.635;      % electrode radius (mm)
Eth = 0.2;        % activation threshold (V/mm)
sig = 0.2;        % isotropic conductivity (S/m)
kap = 4.4;        % anisotropy ratio
sx  = sig*sqrt(kap);     % conductivity along the fibres (x)
sy  = sig/sqrt(kap);     % across the fibres (y).  In 2D the constraint is sx*sy = sig^2
R_ver = 12;                          % domain radius for the verification runs (mm)
hs = [0.4 0.2 0.1 0.05 0.025];       % grid spacings for the convergence studies

% Boundary values: V0 on the electrode (near the centre), 0 on the outer circle
bc = @(x, y) V0 * (x.^2 + y.^2 < 1);

% ----- A. Method of manufactured solutions -----
% V = ln(x^2/sx + y^2/sy) satisfies sx*V_xx + sy*V_yy = 0 exactly,
% so we put it on the boundary and compare the numerical answer with it.
fprintf('\n=== A. Manufactured solution (kappa = %.1f) ===\n', kap);
Vexact = @(x, y) log(x.^2/sx + y.^2/sy);
errM = zeros(size(hs));
for k = 1:length(hs)
    [V, ~, ~, ~, G] = fd_solve3(hs(k), R_ver, a, sx, sy, Vexact);
    errM(k) = max(abs(V(G.dom) - Vexact(G.X(G.dom), G.Y(G.dom))));
end
pM = log2(errM(1:end-1) ./ errM(2:end));     % observed order (h halves each time)
fprintf('h = %6.3f  max error = %.3e\n', hs(1), errM(1));
for k = 2:length(hs)
    fprintf('h = %6.3f  max error = %.3e   observed order p = %.2f\n', hs(k), errM(k), pM(k-1));
end

% ----- B. Isotropic VTA area compared with the analytic answer -----
% |E| = V0/(r ln(R/a))  ->  r_th = V0/(Eth ln(R/a)),  area = pi*r_th^2
fprintf('\n=== B. Isotropic VTA area vs analytic (R = %g mm) ===\n', R_ver);
rth_ex = V0/(Eth*log(R_ver/a));
A_ex = pi*rth_ex^2;                  % includes the electrode (ray method)
A_ex_vox = pi*(rth_ex^2 - a^2);      % excludes the electrode (voxel count)
fprintf('analytic r_th = %.4f mm, area (incl. electrode) = %.4f mm^2\n', rth_ex, A_ex);
Aray = zeros(size(hs)); Avox = Aray; Aray_an = Aray;
for k = 1:length(hs)
    [~, ~, ~, Em, G] = fd_solve3(hs(k), R_ver, a, sig, sig, bc);    % isotropic
    [Aray(k), Avox(k)] = vta_area_2d4(Em, G, Eth);
    [~, ~, ~, Em, G] = fd_solve3(hs(k), R_ver, a, sx, sy, bc);      % anisotropic
    Aray_an(k) = vta_area_2d4(Em, G, Eth);
end
eRay = abs(Aray - A_ex);
eVox = abs(Avox - A_ex_vox);
for k = 1:length(hs)
    fprintf('h = %6.3f  A_ray = %.4f (err %.2e)   A_vox = %.4f (err %.2e)   A_ray(aniso) = %.4f\n', ...
        hs(k), Aray(k), eRay(k), Avox(k), eVox(k), Aray_an(k));
end
% observed order from three successive grids:  p = log2( (A_2h - A_4h) / (A_h - A_2h) )
dA = abs(diff(Aray_an));
p_obs = log2(dA(1:end-1) ./ dA(2:end));
fprintf('anisotropic area, observed orders p from successive triples: %s\n', mat2str(round(p_obs*100)/100));
% Richardson extrapolation with p = 2 (the order confirmed in A)
p = 2;
A_rich = Aray_an(end) + (Aray_an(end) - Aray_an(end-1))/(2^p - 1);
fprintf('Richardson-extrapolated anisotropic area = %.4f mm^2 ; error estimate of finest grid = %.2e\n', ...
    A_rich, abs(A_rich - Aray_an(end)));
c1 = polyfit(log(hs), log(eRay), 1);     % slope of the log-log error plot = order
c2 = polyfit(log(hs), log(eVox), 1);
fprintf('least-squares slope of log(error) vs log(h): ray method %.2f, voxel count %.2f\n', c1(1), c2(1));

% ----- C. Trapezium rule in the angle -----
fprintf('\n=== C. Composite trapezium rule over rays (h = 0.05) ===\n');
[~, ~, ~, Em, G] = fd_solve3(0.05, R_ver, a, sx, sy, bc);
Mlist = [4 8 16 32 64 128];
Aref = vta_area_2d4(Em, G, Eth, 1024);        % reference answer: 1024 rays
eTrap = zeros(size(Mlist));
for k = 1:length(Mlist)
    eTrap(k) = abs(vta_area_2d4(Em, G, Eth, Mlist(k)) - Aref);
    fprintf('M = %4d rays: |A_M - A_1024| = %.2e\n', Mlist(k), eTrap(k));
end

% ----- D. Truncation vs round-off in the central difference -----
% Use the exact V(r), so the only error comes from the difference formula.
fprintf('\n=== D. Central difference for E: truncation vs round-off ===\n');
Rd = 48; r0 = 3;
Vr  = @(r) V0*log(Rd./r)/log(Rd/a);
Eex = V0/(r0*log(Rd/a));
hd  = 10.^(-12:0.25:0);
eD  = abs(-(Vr(r0+hd) - Vr(r0-hd))./(2*hd) - Eex);
M3  = 2*V0/(log(Rd/a)*r0^3);                 % |V'''(r0)|
h_opt = (3*eps*Vr(r0)/M3)^(1/3);
[emin, imin] = min(eD);
fprintf('optimal h (theory) = %.2e mm ; observed minimum error %.2e at h = %.2e mm\n', h_opt, emin, hd(imin));
fprintf('error at the grid spacing used (h = 0.1 mm): %.2e V/mm (truncation dominated)\n', ...
    abs(-(Vr(r0+0.1) - Vr(r0-0.1))/0.2 - Eex));

% ----- E. Condition number of the linear system -----
fprintf('\n=== E. Condition number of the FD matrix (R = %g mm) ===\n', R_ver);
for h = [0.8 0.4 0.2]
    for kk = [1 kap]
        [~, ~, ~, ~, ~, A] = fd_solve3(h, R_ver, a, sig*sqrt(kk), sig/sqrt(kk), @(x,y) 0*x);
        c = condest(A);
        fprintf('h = %.2f, kappa = %.1f : cond_1(A) ~ %.2e ; bound cond*eps = %.1e\n', h, kk, c, c*eps);
    end
end

% ----- F. Domain size -----
fprintf('\n=== F. Domain size (h = 0.2 mm) ===\n');
Rlist = [12 24 48 96];
Ad = zeros(length(Rlist), 2);                % column 1 isotropic, column 2 anisotropic
for k = 1:length(Rlist)
    [~, ~, ~, Em, G] = fd_solve3(0.2, Rlist(k), a, sig, sig, bc);
    Ad(k,1) = vta_area_2d4(Em, G, Eth) - pi*a^2;
    [~, ~, ~, Em, G] = fd_solve3(0.2, Rlist(k), a, sx, sy, bc);
    Ad(k,2) = vta_area_2d4(Em, G, Eth) - pi*a^2;
    A_analytic = pi*((V0/(Eth*log(Rlist(k)/a)))^2 - a^2);
    fprintf('R = %3d mm: iso %.2f mm^2 (analytic %.2f), aniso %.2f mm^2, change %+.1f%%\n', ...
        Rlist(k), Ad(k,1), A_analytic, Ad(k,2), 100*(Ad(k,2)/Ad(k,1) - 1));
end

% ----- G. Main 2D comparison (R = 48 mm, h = 0.1 mm) -----
fprintf('\n=== G. Main 2D comparison (R = 48 mm, h = 0.1 mm, V0 = %g V, kappa = %.1f) ===\n', V0, kap);
R = 48; h = 0.1; Mr = 256;
[~, Exi, Eyi, Emi, G] = fd_solve3(h, R, a, sig, sig, bc);    % isotropic
[~, Exa, Eya, Ema]    = fd_solve3(h, R, a, sx, sy, bc);      % anisotropic
[Ai, ~, ri] = vta_area_2d4(Emi, G, Eth, Mr);   Ai = Ai - pi*a^2;   % area without the electrode
[Aa, ~, ra] = vta_area_2d4(Ema, G, Eth, Mr);   Aa = Aa - pi*a^2;
% spread along x = r_th(0) + r_th(pi);  along y = r_th(pi/2) + r_th(3pi/2)
si = [ri(1) + ri(Mr/2+1), ri(Mr/4+1) + ri(3*Mr/4+1)];
sa = [ra(1) + ra(Mr/2+1), ra(Mr/4+1) + ra(3*Mr/4+1)];
% Dice coefficient between the two activated regions
act_i = (Emi >= Eth) & G.dom;
act_a = (Ema >= Eth) & G.dom;
dice2 = 2*nnz(act_i & act_a) / (nnz(act_i) + nnz(act_a));
% relative field error eps_E, in the ring 1 mm < r < 10 mm
Om = (G.r > 1) & (G.r < 10);
epsE = sqrt(sum((Exa(Om)-Exi(Om)).^2 + (Eya(Om)-Eyi(Om)).^2)) / sqrt(sum(Exi(Om).^2 + Eyi(Om).^2));
fprintf('area: iso %.2f, aniso %.2f mm^2 (%+.1f%%)\n', Ai, Aa, 100*(Aa/Ai-1));
fprintf('spread along fibres (x): iso %.2f, aniso %.2f mm (%+.1f%%)\n', si(1), sa(1), 100*(sa(1)/si(1)-1));
fprintf('spread across fibres (y): iso %.2f, aniso %.2f mm (%+.1f%%)\n', si(2), sa(2), 100*(sa(2)/si(2)-1));
fprintf('aspect ratio x/y: iso %.2f, aniso %.2f ; Dice = %.3f ; eps_E = %.3f\n', si(1)/si(2), sa(1)/sa(2), dice2, epsE);
% aspect ratio as the voltage grows (linear problem: |E| scales with V0)
V0s = [2 3 5 8 12];
asp = zeros(size(V0s));
for k = 1:length(V0s)
    [~, ~, rr] = vta_area_2d4(Ema*V0s(k)/V0, G, Eth, Mr);
    asp(k) = (rr(1) + rr(Mr/2+1)) / (rr(Mr/4+1) + rr(3*Mr/4+1));
end
fprintf('aspect ratio vs V0 %s V: %s  (2D far-field limit of the |E| contour = 1; isopotentials -> sqrt(kappa) = %.2f)\n', ...
    mat2str(V0s), mat2str(round(asp*100)/100), sqrt(kap));
diary off;

% ----- Save results -----
save('Finalresults/results2d.mat', 'hs', 'errM', 'eRay', 'eVox', 'Mlist', 'eTrap', 'hd', 'eD', 'h_opt', ...
     'Aray', 'Avox', 'Aray_an', 'A_rich', 'Rlist', 'Ad', 'Ai', 'Aa', 'si', 'sa', 'dice2', 'epsE', 'V0s', 'asp');

% ----- Figure 3: |E| around the electrode, isotropic vs anisotropic -----
win = abs(G.x) <= 6;                          
x = G.x(win);
Ei = Emi(win, win);
Ea = Ema(win, win);
relE = sqrt((Exa-Exi).^2 + (Eya-Eyi).^2) ./ max(Emi, eps);   
relE(~G.dom) = NaN;
relE = relE(win, win);

figure('Position', [50 50 1200 380]);
subplot(1,3,1);
imagesc(x, x, min(Ei,1)'); axis xy equal tight; hold on;
contour(x, x, Ei', [Eth Eth], 'k', 'LineWidth', 2);
colorbar; caxis([0 1]);
xlabel('x (mm)'); ylabel('y (mm)'); title('(a) Isotropic |E| (V/mm)');

subplot(1,3,2);
imagesc(x, x, min(Ea,1)'); axis xy equal tight; hold on;
contour(x, x, Ea', [Eth Eth], 'k', 'LineWidth', 2);
contour(x, x, Ei', [Eth Eth], 'r--', 'LineWidth', 1.2);
colorbar; caxis([0 1]);
xlabel('x (mm)'); ylabel('y (mm)'); title('(b) Anisotropic, \kappa = 4.4, fibres along x');

subplot(1,3,3);
imagesc(x, x, min(relE,1)'); axis xy equal tight; hold on;
contour(x, x, Ea', [Eth Eth], 'k', 'LineWidth', 1.5);
contour(x, x, Ei', [Eth Eth], 'r--', 'LineWidth', 1.2);
colorbar; caxis([0 1]);
xlabel('x (mm)'); ylabel('y (mm)'); title('(c) |E_{an} - E_{iso}| / |E_{iso}|');
colormap(flipud(bone));
print('Finalresults/fig3_2d_fields.png', '-dpng', '-r200');

disp('part1_2d finished');
