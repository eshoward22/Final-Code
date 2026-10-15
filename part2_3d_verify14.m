% ----- Checking the 3D finite element solver -----

%  (i)   Manufactured solution (isotropic and rotated anisotropic tensor)
%        -> observed order of convergence
%  (ii)  Coaxial cylinders, where the VTA is known exactly
%  (iii) Mesh convergence of the lead-model VTA (observed order, Richardson)
%  (iv)  Domain size: does moving the grounded boundary change the answer?


clear; close all;
if ~exist('Finalresults', 'dir'), mkdir('Finalresults'); end
if exist('Finalresults/summary3d_verify.txt', 'file'), delete('Finalresults/summary3d_verify.txt'); end
diary('Finalresults/summary3d_verify.txt');

Eth = 0.2;
rho = sqrt(2);                       % mesh refinement ratio
hcs = 0.5 ./ rho.^(0:3);             % core spacings 0.5, 0.354, 0.25, 0.177 mm
region = [-6 6 -6 6 -3 7];              % region where the VTA is measured (mm)


% (zcore = [-3 7], Rcore = 6.5, zpad = 30) so the many solves are faster.

% ----- Manufactured solution -----
% With y = S^(-1/2) x, the function u = exp(alpha*y1) cos(alpha*y2) satisfies
% div(S grad u) = 0.  We fix u on the outer boundary and the whole lead
% surface, solve, and compare with u at every node.
fprintf('\n=== (i) Manufactured solution on the lead mesh ===\n');
cases = [1 0; 4.4 45];               % [kappa theta]
errMMS = zeros(length(hcs), 2);
alpha = 0.3;
for k = 1:length(hcs)
    P = dbs_params1(hcs(k));
    P.zcore = [-3 7];  P.Rcore = 6.5;  P.zpad = 30;     % smaller test domain
    P.Rout = P.Rcore;  P.zpad = 0;   % just the core cylinder, no stretched layers
    M = make_mesh6(P);
    for c = 1:2
        S = sigma_tensor2(0.2, cases(c,1), cases(c,2));
        T = inv(sqrtm(S/0.2));       % S^(-1/2) (scaled)
        y1 = T(1,1)*M.X + T(1,2)*M.Y + T(1,3)*(M.Z - 2);
        y2 = T(2,1)*M.X + T(2,2)*M.Y + T(2,3)*(M.Z - 2);
        u = exp(alpha*y1) .* cos(alpha*y2);
        K = assemble_K8(M, S);                       % no encapsulation layer here
        fixed = unique([M.boundary; M.leadsurf]);
        V = solve_fem9(K, fixed, u(fixed), 1e-12);
        errMMS(k,c) = max(abs(V - u)) / max(abs(u));
    end
end
for c = 1:2
    p = log(errMMS(1:end-1,c) ./ errMMS(2:end,c)) / log(rho);
    fprintf('kappa = %.1f, theta = %2d deg: rel. max error %s ; observed p %s\n', cases(c,1), cases(c,2), ...
        mat2str(errMMS(:,c)', 3), mat2str(round(p'*100)/100));
end

% ----- Coaxial cylinders -----
% V = V0 ln(R/r)/ln(R/a),  |E| = V0/(r ln(R/a))
fprintf('\n=== (ii) Coaxial-cylinder analytic test (isotropic, V0 = 2 V) ===\n');
V0 = 2;
coax = zeros(length(hcs), 2);        % columns: VTA cross-section area, r_th
for k = 1:length(hcs)
    P = dbs_params1(hcs(k));
    P.zcore = [-3 7];  P.Rcore = 6.5;  P.zpad = 30;     % smaller test domain
    M = make_mesh6(P);
    rN = sqrt(M.X.^2 + M.Y.^2);
    u = V0*log(P.Rout./rN)/log(P.Rout/P.a);
    K = assemble_K8(M, 0.2*eye(3));
    fixed = unique([M.boundary; M.leadsurf]);
    V = solve_fem9(K, fixed, u(fixed));
    [E, G] = field_on_grid10(M, V, 0.05, [-4 4 -4 4 1 3]);
    area = vta_volume11(E.Em, G.h, Eth) / (G.z(end) - G.z(1) + G.h);   % volume / height
    coax(k,:) = [area, sqrt(area/pi + P.a^2)];
end
rth_ex = V0/(Eth*log(P.Rout/P.a));
fprintf('analytic r_th = %.4f mm\n', rth_ex);
for k = 1:length(hcs)
    fprintf('hc = %.4f: r_th = %.4f mm (error %.2e mm)\n', hcs(k), coax(k,2), abs(coax(k,2) - rth_ex));
end

% ----- Mesh convergence of the lead model -----
% contact c1, 3 V, kappa = 4.4, theta = 90
fprintf('\n=== (iii) Mesh convergence, lead model, contact 1, V0 = 3 V ===\n');
V0 = 3;
vols = zeros(length(hcs), 2);        % columns: isotropic, anisotropic
tensors = {sigma_tensor2(0.2, 1, 0), sigma_tensor2(0.2, 4.4, 90)};
for k = 1:length(hcs)
    P = dbs_params1(hcs(k));
    P.zcore = [-3 7];  P.Rcore = 6.5;  P.zpad = 30;     % smaller test domain
    M = make_mesh6(P);
    fixed = [M.contacts{2}; M.boundary];                          % contact c1 is number 2 in MATLAB
    values = [ones(length(M.contacts{2}),1); zeros(length(M.boundary),1)];   % 1 V and 0 V
    for m = 1:2
        K = assemble_K8(M, tensors{m}, P.sigma_enc);
        [V, info] = solve_fem9(K, fixed, values);
        E = field_on_grid10(M, V, 0.1, region);
        vols(k,m) = vta_volume11(V0*E.Em, 0.1, Eth);     % linear: field at 3 V = 3 * field at 1 V
    end
    fprintf('hc = %.4f mm (%7d unknowns, PCG %3d its, %5.1f s): V_iso = %.2f, V_aniso = %.2f mm^3, change %+.2f%%\n', ...
        hcs(k), info.nfree, info.iter, info.time, vols(k,1), vols(k,2), 100*(vols(k,2)/vols(k,1) - 1));
end
% observed order from successive triples of meshes
pobs = zeros(2, length(hcs)-2);
for m = 1:2
    d = diff(vols(:,m));
    pobs(m,:) = log(abs(d(1:end-1) ./ d(2:end))) / log(rho);
end
fprintf('observed orders p (iso) %s, (aniso) %s\n', mat2str(round(pobs(1,:)*100)/100), mat2str(round(pobs(2,:)*100)/100));
% Richardson extrapolation, using the order from the finest triple (kept between 1 and 2)
pr = min(max(pobs(:,end)', 1), 2);
vr = vols(end,:) + (vols(end,:) - vols(end-1,:)) ./ (rho.^pr - 1);
chg = 100*(vols(:,2)./vols(:,1) - 1);
fprintf('Richardson: V_iso* = %.2f (p = %.2f), V_aniso* = %.2f (p = %.2f), change* = %+.2f%%\n', vr(1), pr(1), vr(2), pr(2), 100*(vr(2)/vr(1) - 1));
fprintf('estimated error of the hc = 0.25 mm mesh: iso %.2f%%, aniso %.2f%%, relative change %.2f percentage points\n', ...
    100*abs(vols(3,1)/vr(1) - 1), 100*abs(vols(3,2)/vr(2) - 1), abs(chg(3) - 100*(vr(2)/vr(1) - 1)));


% ----- Domain size -----
fprintf('\n=== (iv) Domain size (hc = 0.25 mm) ===\n');
Rl = [20 40 80];
vd = zeros(length(Rl), 2);
for k = 1:length(Rl)
    P = dbs_params1(0.25);
    P.zcore = [-3 7];  P.Rcore = 6.5;  P.zpad = 30;     % smaller test domain
    P.Rout = Rl(k);  P.zpad = Rl(k) - 5;
    M = make_mesh6(P);
    fixed = [M.contacts{2}; M.boundary];
    values = [ones(length(M.contacts{2}),1); zeros(length(M.boundary),1)];
    for m = 1:2
        K = assemble_K8(M, tensors{m}, P.sigma_enc);
        V = solve_fem9(K, fixed, values);
        E = field_on_grid10(M, V, 0.1, region);
        vd(k,m) = vta_volume11(3*E.Em, 0.1, Eth);
    end
    fprintf('outer radius %2d mm: V_iso = %.2f, V_aniso = %.2f mm^3, change %+.2f%%\n', Rl(k), vd(k,1), vd(k,2), 100*(vd(k,2)/vd(k,1) - 1));
end
diary off;
save('Finalresults/verify3d.mat', 'hcs', 'errMMS', 'coax', 'rth_ex', 'vols', 'vr', 'pr', 'pobs', 'Rl', 'vd');

% ----- Figure 2 -----
% all the convergence plots (2D results from part 1, 3D from above)
R2 = load('Finalresults/results2d.mat');
figure('Position', [50 50 1000 820]);

subplot(2,2,1);
loglog(R2.hs, R2.errM, 'o-', R2.hs, R2.eRay, 's-', R2.hs, R2.eVox, '^-', 'LineWidth', 1.5); hold on;
loglog(R2.hs, 0.5*R2.errM(1)*(R2.hs/R2.hs(1)).^2, 'k--');
loglog(R2.hs, 0.5*R2.eVox(1)*(R2.hs/R2.hs(1)), 'k:');
legend('MMS potential (max)', 'VTA area: rays + trapezium', 'VTA area: voxel count', 'O(h^2)', 'O(h)', 'Location', 'northwest');
xlabel('grid spacing h (mm)'); ylabel('error'); title('(a) 2D grid convergence'); grid on; ylim([1e-5 1e2]);

subplot(2,2,2);
semilogy(R2.Mlist, R2.eTrap, 'o-', 'LineWidth', 1.5);
xlabel('number of rays M'); ylabel('|A_M - A_{1024}| (mm^2)'); title('(b) Trapezium rule in the angle'); grid on;

subplot(2,2,3);
Vr3 = 2*log(48/3)/log(48/0.635);          % |V| at r0 = 3 mm (as in part 1, D)
M3 = 2*2/(log(48/0.635)*3^3);
loglog(R2.hd, R2.eD, '.-'); hold on;
loglog(R2.hd, M3/6*R2.hd.^2, 'k--');
loglog(R2.hd, eps*Vr3./R2.hd, 'k:');
loglog(R2.h_opt*[1 1], [1e-14 1], 'r-.');
legend('|error|', 'truncation M_3h^2/6', 'round-off \epsilon_{mach}|V|/h', 'h_{opt}', 'Location', 'north');
xlabel('step size h (mm)'); ylabel('|error in E| (V/mm)'); title('(c) Truncation vs round-off'); grid on; ylim([1e-14 1]);

subplot(2,2,4);
loglog(hcs, errMMS(:,1), 'o-', hcs, errMMS(:,2), 's--', hcs, abs(coax(:,2) - rth_ex)/rth_ex, '^-', 'LineWidth', 1.5); hold on;
loglog(hcs, 0.8*errMMS(1,1)*(hcs/hcs(1)).^2, 'k--');
legend('MMS, isotropic', 'MMS, \kappa = 4.4, \theta = 45', 'coaxial r_{th}', 'O(h^2)', 'Location', 'northwest');
xlabel('core mesh spacing h_c (mm)'); ylabel('relative error'); title('(d) 3D finite elements'); grid on;
print('Finalresults/fig2_verification.png', '-dpng', '-r200');

disp('part2_3d_verify finished');
