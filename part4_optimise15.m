% ----- Optimising -----

%  Target T (STN):               ellipsoid, centre (0, 0.75, 3) mm, semi-axes (2.5, 2.0, 3.0) mm
%  Avoidance A (internal capsule): slab  -4 <= y <= -2.5,  |x| <= 4,  -1 <= z <= 7 mm
%  Choices: contact c in {0,1,2,3} and voltage V0 = 0.5, 0.55, ..., 5 V (91 values)
%  Objective:  J = coverage - lambda*spill
%       coverage = |VTA n T| / |T|,     spill = |VTA n A| / |A|
%  Because the problem is linear, VTA(c, V0) = { |E1(c)| >= Eth/V0 } where E1 is
%  the field for 1 V, so we only need ONE solve per contact and tissue model.

clear; close all;
if ~exist('Finalresults', 'dir'), mkdir('Finalresults'); end
if exist('Finalresults/summary3d_opt.txt', 'file'), delete('Finalresults/summary3d_opt.txt'); end
diary('Finalresults/summary3d_opt.txt');

P = dbs_params1(0.25);
M = make_mesh6(P);
Eth = P.Eth;  h = 0.1;  region = [-6 6 -6 6 -4 10];
names  = {'isotropic', 'kappa 4.4, theta 0', 'kappa 4.4, theta 45', 'kappa 4.4, theta 90', 'kappa 10, theta 90'};
kappas = [1 4.4 4.4 4.4 10];
thetas = [0 0   45  90  90];
V0grid  = 0.5:0.05:5;
lambdas = [0.5 1 2];               % lambda = 1 is the main setting, 0.5 and 2 test sensitivity
nm = 5;  nc = 4;  nv = length(V0grid);

% ----- Step 1 -----
% Coverage and spill for every model, contact and voltage
% arrays are indexed (model, contact, voltage)
cov = zeros(nm, nc, nv);  spill = cov;  pout = cov;  vvol = cov;
for m = 1:nm
    K = assemble_K8(M, sigma_tensor2(P.sigma_iso, kappas(m), thetas(m)), P.sigma_enc);
    for c = 1:nc
        fixed  = [M.contacts{c}; M.boundary];                       % contact c at 1 V, boundary at 0 V
        values = [ones(length(M.contacts{c}),1); zeros(length(M.boundary),1)];
        V = solve_fem9(K, fixed, values);
        [E, G] = field_on_grid10(M, V, h, region);

        if m == 1 && c == 1           % build the regions T and A (only needed once)
            [X, Y, Z] = ndgrid(G.x, G.y, G.z);
            T = ((X/2.5).^2 + ((Y-0.75)/2.0).^2 + ((Z-3)/3.0).^2 <= 1) & ~G.lead;
            A = (Y >= -4) & (Y <= -2.5) & (abs(X) <= 4) & (Z >= -1) & (Z <= 7);
            tissue = ~G.lead;
            fprintf('|T| = %.1f mm^3, |A| = %.1f mm^3 (voxel counts)\n', nnz(T)*h^3, nnz(A)*h^3);
        end

        eT = E.Em(T);  eA = E.Em(A);  eAll = E.Em(tissue);   % the 1 V field in T, in A, everywhere
        for k = 1:nv
            thr = Eth / V0grid(k);                % VTA at V0  =  { |E1| >= Eth/V0 }
            nT = nnz(eT >= thr);                  % activated voxels in T
            nA = nnz(eA >= thr);                  % activated voxels in A
            nV = nnz(eAll >= thr);                % activated voxels in total
            cov(m,c,k)   = nT / length(eT);
            spill(m,c,k) = nA / length(eA);
            vvol(m,c,k)  = nV * h^3;
            pout(m,c,k)  = 1 - nT/max(nV, 1);     % fraction of the VTA outside T
        end
    end
    fprintf('model %-20s solved\n', names{m});
end

% ----- Step 2 -----
% Exhaustive search over the 4 x 91 = 364 settings
best  = zeros(nm, 3, 2);           % (model, lambda, [contact, voltage index])
Jbest = zeros(nm, 3);
for m = 1:nm
    for l = 1:3
        Jbest(m,l) = -Inf;
        for k = 1:nv               % lowest voltage first, so a tie goes to the lower voltage
            for c = 1:nc
                J = cov(m,c,k) - lambdas(l)*spill(m,c,k);
                if J > Jbest(m,l)
                    Jbest(m,l) = J;
                    best(m,l,:) = [c, k];
                end
            end
        end
    end
end
fprintf('\n=== Optimal settings (lambda = 1) ===\n');
l0 = 2;                            % lambdas(2) = 1
for m = 1:nm
    c = best(m,l0,1);  k = best(m,l0,2);
    fprintf('%-20s: contact %d, V0 = %.2f V, coverage %.1f%%, spill %.1f%%, overspill(Cubo) %.1f%%, VTA %.1f mm^3, J = %.3f\n', ...
        names{m}, c-1, V0grid(k), 100*cov(m,c,k), 100*spill(m,c,k), 100*pout(m,c,k), vvol(m,c,k), Jbest(m,l0));
end

% ----- Step 3 -----
% Model mismatch - use the isotropic optimum in anisotropic tissue
fprintf('\n=== Isotropic optimum evaluated in anisotropic tissue ===\n');
mis = zeros(nm, 3, 4);             % [coverage, spill, J, loss compared with the anisotropic optimum]
for l = 1:3
    ci = best(1,l,1);  ki = best(1,l,2);          % isotropic optimum
    for m = 2:nm
        Jm = cov(m,ci,ki) - lambdas(l)*spill(m,ci,ki);
        mis(m,l,:) = [cov(m,ci,ki), spill(m,ci,ki), Jm, Jbest(m,l) - Jm];
        fprintf('lambda %.1f, %-20s: iso optimum (c%d, %.2f V) gives coverage %.1f%%, spill %.1f%%, J = %.3f; anisotropic optimum J = %.3f (loss %.3f)\n', ...
            lambdas(l), names{m}, ci-1, V0grid(ki), 100*mis(m,l,1), 100*mis(m,l,2), mis(m,l,3), Jbest(m,l), mis(m,l,4));
    end
end
% the voltage at which spill first goes above 1% on the isotropic-optimal contact
ci = best(1,l0,1);
for m = 1:nm
    k1 = find(squeeze(spill(m,ci,:)) > 0.01, 1);
    if isempty(k1)
        fprintf('%-20s: spill never exceeds 1%% on contact %d\n', names{m}, ci-1);
    else
        fprintf('%-20s: spill exceeds 1%% on contact %d above V0 = %.2f V\n', names{m}, ci-1, V0grid(k1));
    end
end
diary off;
save('Finalresults/opt3d.mat', 'names', 'V0grid', 'lambdas', 'cov', 'spill', 'pout', 'vvol', 'best', 'Jbest', 'mis');


% ----- Figure 6 -----
% Coverage vs spill, and the objective, on the isotropic-optimal contact
ci = best(1,l0,1);  ki = best(1,l0,2);           % isotropic optimum
ma = 4;                                          % anisotropic model kappa 4.4, theta 90
ca = best(ma,l0,1); ka = best(ma,l0,2);          % its optimum
figure('Position', [50 50 1000 400]);
subplot(1,2,1); hold on;
plot(100*squeeze(spill(1,ci,:)),  100*squeeze(cov(1,ci,:)),  'b-', 'LineWidth', 1.8);
plot(100*squeeze(spill(ma,ci,:)), 100*squeeze(cov(ma,ci,:)), 'r-', 'LineWidth', 1.8);
plot(100*spill(1,ci,ki),  100*cov(1,ci,ki),  'bo', 'MarkerFaceColor', 'b', 'MarkerSize', 8);
plot(100*spill(ma,ca,ka), 100*cov(ma,ca,ka), 'rs', 'MarkerFaceColor', 'r', 'MarkerSize', 8);
plot(100*spill(ma,ci,ki), 100*cov(ma,ci,ki), 'kx', 'MarkerSize', 12, 'LineWidth', 2);
legend(sprintf('isotropic, contact %d', ci-1), sprintf('anisotropic, contact %d', ci-1), ...
       'isotropic optimum', 'anisotropic optimum', 'isotropic optimum in anisotropic tissue', 'Location', 'southeast');
xlabel('spill into avoidance region A (%)'); ylabel('coverage of target T (%)'); grid on; box on;
title('(a) Coverage vs spill as V_0 increases');
subplot(1,2,2); hold on;
Ji = squeeze(cov(1,ci,:)  - spill(1,ci,:));
Ja = squeeze(cov(ma,ci,:) - spill(ma,ci,:));
plot(V0grid, Ji, 'b-', V0grid, Ja, 'r-', 'LineWidth', 1.8);
plot(V0grid(ki), Ji(ki), 'bo', 'MarkerFaceColor', 'b', 'MarkerSize', 8);
plot(V0grid(ki), Ja(ki), 'kx', 'MarkerSize', 12, 'LineWidth', 2);
legend('isotropic model', 'anisotropic model (\kappa = 4.4, \theta = 90)', 'Location', 'southwest');
xlabel('V_0 (V)'); ylabel('J = coverage - spill'); grid on; box on;
title(sprintf('(b) Objective on contact %d', ci-1));
print('Finalresults/fig6_optimisation.png', '-dpng', '-r200');


% ----- Figure 1 -----
% The model geometry (lead, contacts, T and A), to scale
t = linspace(0, 2*pi, 200);
figure('Position', [50 50 900 430]);
subplot(1,2,1); hold on;                                    % plane x = 0
fill([-4 -2.5 -2.5 -4], [-1 -1 7 7], [1 0.85 0.85]);        % A
fill([-P.a P.a P.a -P.a], [-4 -4 10 10], [0.55 0.55 0.55]); % lead
for c = 1:4
    fill([-P.a P.a P.a -P.a], P.zc(c) + [-0.75 -0.75 0.75 0.75], [0.15 0.15 0.15]);
    text(-P.a - 1.7, P.zc(c), sprintf('c%d', c-1));
end
plot(0.75 + 2.0*cos(t), 3 + 3.0*sin(t), 'g-', 'LineWidth', 2);   % T
text(-3.9, 7.6, 'A'); text(1.9, 6.3, 'T');
axis equal; xlim([-6 6]); ylim([-4.5 9]); grid on;
xlabel('y (mm)'); ylabel('z (mm)'); title('(a) Plane x = 0 (lead along z)');
subplot(1,2,2); hold on;                                    % plane z = 3
fill([-4 4 4 -4], [-4 -4 -2.5 -2.5], [1 0.85 0.85]);         % A
plot(2.5*cos(t), 0.75 + 2.0*sin(t), 'g-', 'LineWidth', 2);    % T
fill((P.a + P.t_enc)*cos(t), (P.a + P.t_enc)*sin(t), [0.85 0.85 0.85]);   % encapsulation
fill(P.a*cos(t), P.a*sin(t), [0.15 0.15 0.15]);                           % lead
text(-3.8, -3.3, 'A'); text(1.9, 2.6, 'T');
axis equal; xlim([-6 6]); ylim([-6 6]); grid on;
xlabel('x (mm)'); ylabel('y (mm)'); title('(b) Axial plane z = 3 mm');
print('Finalresults/fig1_geometry.png', '-dpng', '-r200');

disp('part4_optimise finished');
