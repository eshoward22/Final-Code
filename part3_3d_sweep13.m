% ----- Isotropic vs anisotropic  -----

%  Contact c1 (centre z = 2 mm), V0 = 3 V, Eth = 0.2 V/mm.
%  kappa in {1, 2, 4.4, 7, 10}, theta in {0, 45, 90} degrees.
%  For every (kappa, theta): VTA volume, spreads, sphericity, Dice, field error eps_E.
%  Then a spline + bisection finds where the H1 thresholds are crossed.

clear; close all;
if ~exist('Finalresults', 'dir'), mkdir('Finalresults'); end
if exist('Finalresults/summary3d_sweep.txt', 'file'), delete('Finalresults/summary3d_sweep.txt'); end
diary('Finalresults/summary3d_sweep.txt');

P = dbs_params1(0.25);
M = make_mesh6(P);
Eth = P.Eth;  V0 = 3;
region = [-6 6 -6 6 -3 7];  h = 0.1;          % where the VTA is measured, and the grid spacing (mm)
fprintf('mesh: %d nodes, %d elements\n', length(M.X), size(M.conn,1));

% Boundary conditions: contact c1 (number 2 in MATLAB) at 1 V, outer boundary at 0 V.
% The problem is linear, so the field at V0 volts is just V0 * (field at 1 V).
fixed  = [M.contacts{2}; M.boundary];
values = [ones(length(M.contacts{2}),1); zeros(length(M.boundary),1)];

% ----- Isotropic reference -----
K = assemble_K8(M, sigma_tensor2(P.sigma_iso, 1, 0), P.sigma_enc);
[V, info] = solve_fem9(K, fixed, values);
[Ei, G] = field_on_grid10(M, V, h, region);
Ri = vta_shape12(V0*Ei.Em, G, Eth);
[Xs, Ys, ~] = ndgrid(G.x, G.y, G.z);
Om = sqrt(Xs.^2 + Ys.^2) > P.a + P.t_enc + 0.05;   % region for eps_E: outside the encapsulation layer
fprintf('isotropic: V = %.2f mm^3, spreads x/y/z = %.1f/%.1f/%.1f mm, sphericity %.3f (PCG %d its, %.1f s)\n', ...
    Ri.volume, Ri.spread, Ri.sphericity, info.iter, info.time);

% ----- Sweep over kappa and theta -----
kap = [1 2 4.4 7 10];
th  = [0 45 90];
vol = zeros(length(kap), length(th));  dice = vol;  sph = vol;  epsE = vol;
spr = zeros(length(kap), length(th), 3);
% kappa = 1 is isotropic for every theta, so row 1 is just the reference
vol(1,:) = Ri.volume;  dice(1,:) = 1;  sph(1,:) = Ri.sphericity;  epsE(1,:) = 0;
for b = 1:3
    spr(1,b,:) = Ri.spread;
end
for a = 2:length(kap)
    for b = 1:length(th)
        K = assemble_K8(M, sigma_tensor2(P.sigma_iso, kap(a), th(b)), P.sigma_enc);
        [V, info] = solve_fem9(K, fixed, values);
        Ea = field_on_grid10(M, V, h, region);
        Ra = vta_shape12(V0*Ea.Em, G, Eth, V0*Ei.Em);
        vol(a,b) = Ra.volume;  dice(a,b) = Ra.dice;  sph(a,b) = Ra.sphericity;  spr(a,b,:) = Ra.spread;
        % relative field error  eps_E = ||E_an - E_iso|| / ||E_iso||  over Om
        diff2 = (Ea.Ex - Ei.Ex).^2 + (Ea.Ey - Ei.Ey).^2 + (Ea.Ez - Ei.Ez).^2;
        size2 = Ei.Ex.^2 + Ei.Ey.^2 + Ei.Ez.^2;
        epsE(a,b) = sqrt(sum(diff2(Om)) / sum(size2(Om)));
        fprintf('kappa %4.1f theta %2d: V = %6.2f mm^3 (%+5.1f%%), spreads %.1f/%.1f/%.1f mm, Psi %.3f, Dice %.3f, eps_E %.3f, PCG %d its\n', ...
            kap(a), th(b), vol(a,b), 100*(vol(a,b)/Ri.volume - 1), Ra.spread, Ra.sphericity, Ra.dice, epsE(a,b), info.iter);

        if kap(a) == 4.4                             
            if th(b) == 0,  E_t0  = V0*Ea.Em; end
            if th(b) == 45, E_t45 = V0*Ea.Em; end
            if th(b) == 90
                E_t90 = V0*Ea.Em;
             
                V0s = [1 2 3 4 5 6];
                asp = zeros(size(V0s));
                for k = 1:length(V0s)
                    Rk = vta_shape12(V0s(k)*Ea.Em, G, Eth);
                    asp(k) = Rk.spread(2) / Rk.spread(1);
                end
                fprintf('aspect ratio y/x (theta 90, kappa 4.4) for V0 = %s V: %s ; far-field limit kappa^(1/4) = %.2f\n', ...
                    mat2str(V0s), mat2str(round(asp*100)/100), 4.4^(1/4));
            end
        end
    end
end
dV = 100*(vol/Ri.volume - 1);                         % volume change in %


% ----- H1 thresholds -----
% Cubic spline through the 5 kappa values + bisection
% kappa* where the volume change reaches 10%, and where Dice falls to 0.9
kstar = zeros(2, 3);
for b = 1:3
    fV = @(k) interp1(kap, dV(:,b), k, 'spline') - 10;
    fD = @(k) interp1(kap, dice(:,b), k, 'spline') - 0.9;
    kstar(1,b) = bisection(fV, 1, 10);
    kstar(2,b) = bisection(fD, 1, 10);
    fprintf('theta %2d: volume change reaches 10%% at kappa* = %s ; Dice < 0.9 for kappa > %s\n', th(b), ...
        num2str(kstar(1,b), 3), num2str(kstar(2,b), 3));
end
diary off;
save('Finalresults/sweep3d.mat', 'kap', 'th', 'vol', 'dV', 'dice', 'sph', 'epsE', 'spr', 'kstar', 'V0s', 'asp');

% ----- Figure 4 -----
% Outline of the VTA (|E| = Eth) in two planes
E0 = V0*Ei.Em;
iz = find(abs(G.z - 2) < 1e-9);       % plane z = 2 mm (centre of contact c1)
ix = find(abs(G.x) < 1e-9);           % plane x = 0
figure('Position', [50 50 1000 430]);
subplot(1,2,1); hold on;
contour(G.x, G.y, E0(:,:,iz)',    [Eth Eth], 'k',  'LineWidth', 2);
contour(G.x, G.y, E_t0(:,:,iz)',  [Eth Eth], 'b',  'LineWidth', 1.5);
contour(G.x, G.y, E_t45(:,:,iz)', [Eth Eth], 'g',  'LineWidth', 1.5);
contour(G.x, G.y, E_t90(:,:,iz)', [Eth Eth], 'r',  'LineWidth', 1.5);
t = linspace(0, 2*pi, 100);
fill(P.a*cos(t), P.a*sin(t), [0.2 0.2 0.2]);
legend('isotropic', '\kappa = 4.4, \theta = 0', '\kappa = 4.4, \theta = 45', '\kappa = 4.4, \theta = 90');
axis equal; xlim([-5 7.5]); ylim([-5 5]); grid on; box on;
xlabel('x (mm)'); ylabel('y (mm)'); title('(a) Axial plane z = 2 mm');

subplot(1,2,2); hold on;
contour(G.y, G.z, squeeze(E0(ix,:,:))',    [Eth Eth], 'k', 'LineWidth', 2);
contour(G.y, G.z, squeeze(E_t0(ix,:,:))',  [Eth Eth], 'b', 'LineWidth', 1.5);
contour(G.y, G.z, squeeze(E_t45(ix,:,:))', [Eth Eth], 'g', 'LineWidth', 1.5);
contour(G.y, G.z, squeeze(E_t90(ix,:,:))', [Eth Eth], 'r', 'LineWidth', 1.5);
fill([-P.a P.a P.a -P.a], [-3 -3 7 7], [0.6 0.6 0.6]);             % the lead
for c = 1:4
    fill([-P.a P.a P.a -P.a], P.zc(c) + [-0.75 -0.75 0.75 0.75], [0.15 0.15 0.15]);   % contacts
end
axis equal; xlim([-5 5]); ylim([-3 7]); grid on; box on;
xlabel('y (mm)'); ylabel('z (mm)'); title('(b) Plane x = 0 (along the lead)');
print('Finalresults/fig4_crosssections.png', '-dpng', '-r200');

% ----- Figure 5 -----
% Sensitivity to kappa and theta
figure('Position', [50 50 1200 360]);
subplot(1,3,1);
plot(kap, dV, 'o-', 'LineWidth', 1.5); hold on;
plot([1 10], [10 10], 'k--');
legend('\theta = 0', '\theta = 45', '\theta = 90', 'H1: 10%', 'Location', 'northwest');
xlabel('anisotropy ratio \kappa'); ylabel('VTA volume change (%)'); grid on; title('(a) Volume');
subplot(1,3,2);
plot(kap, dice, 'o-', 'LineWidth', 1.5); hold on;
plot([1 10], [0.9 0.9], 'k--');
xlabel('\kappa'); ylabel('Dice coefficient vs isotropic'); grid on; title('(b) Overlap');
subplot(1,3,3);
plot(V0s, asp, 'o-', 'LineWidth', 1.5); hold on;
plot([V0s(1) V0s(end)], 4.4^(1/4)*[1 1], 'k--');
legend('3D model', '\kappa^{1/4} (point-source far field)', 'Location', 'southeast');
xlabel('V_0 (V)'); ylabel('spread along / across fibres'); grid on; title('(c) Elongation, \kappa = 4.4, \theta = 90');
print('Finalresults/fig5_sensitivity.png', '-dpng', '-r200');

disp('part3_3d_sweep finished');
