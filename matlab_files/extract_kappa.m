%% extract_kappa_sweep.m
%
%  Sweeps the edge-to-edge gap between the two waveguide cores,
%  calls compute_modes on the coupled cross-section at each gap,
%  prints a full diagnostic table of ALL modes found, then extracts
%  kappa_TE and kappa_TM and fits K1/K2.
%
%  Run this script FIRST to verify mode ordering before running the
%  full ring_resonator_analysis.m pipeline.

clc; clear; close all;

%% ── User parameters (mirror ring_resonator_analysis.m) ───────────────────
lambda0_um = 0.750;

n_core = get_refractive_index('si3n4', lambda0_um);
n_clad = get_refractive_index('sio2',  lambda0_um);
n_box  = n_clad;

w_core = 0.890;
h_core = 0.670;
h_clad = 2.0;
h_box  = 1.5;
w_sim  = 8.0;

gap_scan_nm = 100:50:800;
gap_scan_um = gap_scan_nm / 1e3;

% Number of modes to request — increase if TM modes are being missed
num_modes_request = 6;

mesh_res = 60;

%% ── Storage ───────────────────────────────────────────────────────────────
k0       = 2*pi / lambda0_um;
N        = numel(gap_scan_um);
kappa_TE = NaN(1, N);
kappa_TM = NaN(1, N);

fprintf('\n%s\n', repmat('─',1,72));
fprintf('  gap(nm)  | modes found |  TE modes (n_eff)         |  TM modes (n_eff)\n');
fprintf('%s\n', repmat('─',1,72));

for ig = 1:N
    gap_um = gap_scan_um(ig);

    [nodes, elems, epsilon_r, ~] = build_coupled_wg_mesh( ...
        w_core, h_core, gap_um, ...
        h_clad, h_box,  w_sim, ...
        n_core, n_clad, n_box, ...
        mesh_res);

    modes = compute_modes(nodes, elems, epsilon_r, lambda0_um, ...
                          'num_modes', num_modes_request, 'mu_r', 1.0);

    neff_all = [modes.n_eff];
    te_frac  = [modes.te_fraction];
    tm_frac  = [modes.tm_fraction];

    % ── Classify modes ────────────────────────────────────────────────────
    % Use a threshold slightly below 0.5 to catch weakly hybrid modes.
    % Adjust TE_THRESH if your waveguide is strongly hybrid.
    TE_THRESH = 0.5;

    te_mask = te_frac >= TE_THRESH;
    tm_mask = ~te_mask;          % everything else treated as TM

    nTE = sort(neff_all(te_mask), 'descend');
    nTM = sort(neff_all(tm_mask), 'descend');

    % ── Diagnostic print ─────────────────────────────────────────────────
    te_str = sprintf('%.5f  ', nTE);
    tm_str = sprintf('%.5f  ', nTM);
    fprintf('  %6.0f   |    %2d       |  %-26s|  %s\n', ...
            gap_scan_nm(ig), numel(modes), te_str, tm_str);

    % ── Extract kappa ─────────────────────────────────────────────────────
    if numel(nTE) >= 2
        kappa_TE(ig) = (nTE(1) - nTE(2)) * k0 / 2;   % [rad/µm]
    end
    if numel(nTM) >= 2
        kappa_TM(ig) = (nTM(1) - nTM(2)) * k0 / 2;
    end
end
fprintf('%s\n', repmat('─',1,72));

%% ── Convert to SI ─────────────────────────────────────────────────────────
gap_m        = gap_scan_um * 1e-6;
kappa_TE_m   = kappa_TE    * 1e6;   % rad/m
kappa_TM_m   = kappa_TM    * 1e6;

fprintf('\nRaw kappa values:\n');
fprintf('  gap(nm)   kappa_TE (rad/m)   kappa_TM (rad/m)\n');
for ig = 1:N
    fprintf('  %6.0f    %+.4e         %+.4e\n', ...
            gap_scan_nm(ig), kappa_TE_m(ig), kappa_TM_m(ig));
end

%% ── Robust exponential fit ────────────────────────────────────────────────
% Only fit where we have valid, positive kappa values (negative kappa
% would indicate a mode-ordering or sign error — flag it).

[fit_TE, K1_TE, K2_TE] = robust_exp_fit(gap_m, kappa_TE_m, 'TE');
[fit_TM, K1_TM, K2_TM] = robust_exp_fit(gap_m, kappa_TM_m, 'TM');

%% ── Plot ──────────────────────────────────────────────────────────────────
d_plot_nm = linspace(gap_scan_nm(1), gap_scan_nm(end), 400);
d_plot_m  = d_plot_nm * 1e-9;

figure('Name','kappa(d) sweep','Color','w','Position',[80 80 720 440]);

ax = gca;
valid_TE = ~isnan(kappa_TE_m) & kappa_TE_m > 0;
valid_TM = ~isnan(kappa_TM_m) & kappa_TM_m > 0;

if any(valid_TE)
    semilogy(gap_scan_nm(valid_TE), kappa_TE_m(valid_TE), 'bo', ...
             'MarkerSize',7,'DisplayName','TE  (FEM)'); hold on;
end
if any(valid_TM)
    semilogy(gap_scan_nm(valid_TM), kappa_TM_m(valid_TM), 'rs', ...
             'MarkerSize',7,'DisplayName','TM  (FEM)'); hold on;
end
if ~isnan(K1_TE)
    semilogy(d_plot_nm, K1_TE*exp(K2_TE*d_plot_m), 'b-', ...
             'LineWidth',1.6,'DisplayName', ...
             sprintf('TE fit  K1=%.2e  K2=%.2e', K1_TE, K2_TE));
end
if ~isnan(K1_TM)
    semilogy(d_plot_nm, K1_TM*exp(K2_TM*d_plot_m), 'r-', ...
             'LineWidth',1.6,'DisplayName', ...
             sprintf('TM fit  K1=%.2e  K2=%.2e', K1_TM, K2_TM));
end

xlabel('Edge-to-edge gap  (nm)');
ylabel('\kappa  (rad m^{-1})');
title('\kappa(d): supermode splitting — FEM coupled cross-section');
legend('Location','northeast'); grid on;

%% ── Save fit coefficients to workspace for main script ───────────────────
assignin('base','K1_TE', K1_TE);
assignin('base','K2_TE', K2_TE);
assignin('base','K1_TM', K1_TM);
assignin('base','K2_TM', K2_TM);
fprintf('\nK1/K2 saved to workspace. Run ring_resonator_analysis.m next.\n');

%% ── Local helper ──────────────────────────────────────────────────────────
function [fo, K1, K2] = robust_exp_fit(gap_m, kappa_m, label)
    valid = ~isnan(kappa_m) & kappa_m > 0;
    n_ok  = sum(valid);

    if n_ok < 2
        fprintf('\n  WARNING (%s): only %d valid kappa points found.\n', label, n_ok);
        fprintf('  Possible causes:\n');
        fprintf('    1. num_modes_request too low → TM supermodes not computed.\n');
        fprintf('    2. TE_THRESH misclassifying TM modes as TE.\n');
        fprintf('    3. Waveguide is strongly TE-like; TM modes are leaky/absent.\n');
        fprintf('  Action: increase num_modes_request (try 8 or 10) and re-run.\n\n');
        fo = []; K1 = NaN; K2 = NaN;
        return
    end

    fo  = fit(gap_m(valid)', kappa_m(valid)', 'exp1');
    K1  = fo.a;
    K2  = fo.b;
    r2  = 1 - sum((kappa_m(valid) - fo(gap_m(valid))').^2) / ...
              sum((kappa_m(valid) - mean(kappa_m(valid))).^2);

    fprintf('\n  %s fit:  K1 = %.4e m⁻¹,  K2 = %.4e m⁻¹,  R² = %.5f\n', ...
            label, K1, K2, r2);

    if K2 >= 0
        fprintf('  WARNING (%s): K2 = %.3e ≥ 0 — kappa INCREASES with gap.\n', label, K2);
        fprintf('  This indicates a mode-ordering error. Check the diagnostic table.\n');
    end
end