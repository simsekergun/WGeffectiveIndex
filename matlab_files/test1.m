%% test 1
% FEM solver for EM modes in dielectric waveguides
% Original formulation is based on: https://doi.org/10.1080/02726340290084012
% here: we apply PML instead of PEC
%
% Solves the full-vectorial eigenvalue problem for guided modes using
% mixed Nedelec (curl-conforming) + Lagrange (nodal) elements on triangles.
%
% Physics: Hybrid-mode formulation with transverse Nedelec and z-component
% Lagrange DOFs, following the weak form from the reference paper.
%
% Usage:
%   Run the script directly for a demo 
%   at 1550 nm, or call compute_modes() from your own script.
%
% Requirements: MATLAB R2020b+ (for complex sparse eigs)

%%
% Resolution 100 ppw
clear; close all;

%% ============================================================
%  TEST 1
% COMSOL Multiphysics
% 1.793112981078253
% 1.7519983832021118
% 1.646924391702016
% 1.6269670822697695
% 
% MATLAB Resolution 100, Elapsed time is 9.291088 seconds.
% Mode 1:  n_eff = 1.792937 + 0.00e+00i,  TE-frac = 0.999
% Mode 2:  n_eff = 1.751915 + 0.00e+00i,  TE-frac = 0.004
% Mode 3:  n_eff = 1.646548 + 0.00e+00i,  TE-frac = 0.993
% Mode 4:  n_eff = 1.626750 + 0.00e+00i,  TE-frac = 0.021
% Mode 5:  n_eff = 1.455926 + 0.00e+00i,  TE-frac = 0.959
% Mode 6:  n_eff = 1.452024 + 0.00e+00i,  TE-frac = 0.899
% 
% MATLAB Resolution 200, Elapsed time is 49.791220 seconds.
% --- Guided modes ---
% Mode 1:  n_eff = 1.793075 + 0.00e+00i,  TE-frac = 0.999
% Mode 2:  n_eff = 1.751982 + 0.00e+00i,  TE-frac = 0.004
% Mode 3:  n_eff = 1.646853 + 0.00e+00i,  TE-frac = 0.993
% Mode 4:  n_eff = 1.626931 + 0.00e+00i,  TE-frac = 0.021
% Mode 5:  n_eff = 1.456183 + 0.00e+00i,  TE-frac = 0.959
% Mode 6:  n_eff = 1.452068 + 0.00e+00i,  TE-frac = 0.900
% 
% MATLAB Resolution 300, Elapsed time is 121.673811 seconds.
% Mode 1:  n_eff = 1.793095 + 0.00e+00i,  TE-frac = 0.999
% Mode 2:  n_eff = 1.751991 + 0.00e+00i,  TE-frac = 0.004
% Mode 3:  n_eff = 1.646890 + 0.00e+00i,  TE-frac = 0.993
% Mode 4:  n_eff = 1.626950 + 0.00e+00i,  TE-frac = 0.021
% Mode 5:  n_eff = 1.456216 + 0.00e+00i,  TE-frac = 0.959
% Mode 6:  n_eff = 1.452073 + 0.00e+00i,  TE-frac = 0.900
%% ============================================================

wavelength = 1.55; % µm
% ---- Geometry (all lengths in micrometres) ----
w_core   = 1.6;   % core width
h_core   = 0.7;   % core height
h_clad   = 2.7;   % cladding height above core
h_box    = 2.00;   % buried-oxide thickness
w_sim    = 6.00;   % total simulation width
n_core   = 1.9761; % refractive index of the core
n_clad   = 1.444;  % refractive index of the cladding
n_box    = n_clad;  % refractive index of the substrate

num_modes  = 6;
mesh_res   = round(w_sim/wavelength*100);   % approximate triangles along the longest side

compute_overlaps = 0; % set it to 1, if you want to calculate the mode overlap matrix
plot_mesh = 0;

fprintf('Building mesh ...\n');
[nodes, elems, epsilon_r, regions] = build_soi_mesh( ...
    w_core, h_core, h_clad, h_box, w_sim, ...
    n_core, n_clad, n_box, mesh_res);

fprintf('Nodes: %d,  Elements: %d\n', size(nodes,1), size(elems,1));


if plot_mesh ==1
    figure('Name', 'Mesh Visualization (Permittivity)', 'Color', 'w');
    ph = patch('Faces', elems, 'Vertices', nodes, ...
               'FaceVertexCData', epsilon_r, ...
               'FaceColor', 'flat', ...
               'EdgeColor', 'k', ...
               'LineWidth', 0.1);
    axis equal tight;
    xlabel('x [\mum]'); ylabel('y [\mum]');
    title('Zoomed in Mesh Material Distribution (\epsilon_r)');
    cb = colorbar;
    ylabel(cb, 'Relative Permittivity \epsilon_r');
    xlim([-w_core*1.2, w_core*1.2]); 
    ylim([-h_core*0.2, h_core*1.2]);
    print -dpng figure_mesh
end

tic
fprintf('Assembling & solving eigenvalue problem ...\n');
modes = compute_modes(nodes, elems, epsilon_r, wavelength, ...
    'num_modes', num_modes, 'mu_r', 1.0);
toc
keyboard
%% ---- Print results ----
fprintf('\n--- Guided modes ---\n');
for m = 1:length(modes)
    fprintf('Mode %d:  n_eff = %.6f + %.2ei,  TE-frac = %.3f\n', ...
        m, real(modes(m).n_eff), imag(modes(m).n_eff), modes(m).te_fraction);
end

%% ---- Plot dominant mode ----
plot_mode_fields(modes(1), nodes, elems, 'Mode 1 (fundamental TE)');
plot_mode_fields(modes(2), nodes, elems, 'Mode 2 ');
plot_mode_fields(modes(3), nodes, elems, 'Mode 3');

if compute_overlaps ==1
%% ---- Overlap matrix ----
    fprintf('\nComputing overlap matrix ...\n');
    N = length(modes);
    OL = zeros(N);
    for i = 1:N
        for j = 1:N
            OL(i,j) = calculate_overlap(modes(i), modes(j));
        end
    end
    figure('Name','Overlap matrix');
    imagesc(real(OL)); colorbar; axis square;
    title('Re\{Overlap integrals\}');
    xlabel('Mode j'); ylabel('Mode i');
    set(gca,'XTick',1:N,'YTick',1:N);
end


