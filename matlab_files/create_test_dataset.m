clear;

%% Parallel setup
parpool_no = 8;

p = gcp('nocreate');
if isempty(p)
    parpool(parpool_no);
elseif p.NumWorkers ~= parpool_no
    delete(gcp('nocreate'));
    parpool(parpool_no);
end

testing_data = zeros(1000, 5);

num_modes  = 2;
compute_overlaps = 0; % set it to 1, if you want to calculate the mode overlap matrix
plot_mesh = 0;
h_clad   = 3;   % cladding height above core
h_box    = 3.00;   % buried-oxide thickness
w_sim    = 6.00;   % total simulation width

parfor ii =1:1000
    wavelength = randi([750, 1650])*1e-3;
    wg_width = randi([750, 3000])*1e-3;
    wg_height = randi([600, 1000])*1e-3;

    n_core = get_refractive_index('si3n4', wavelength);
    n_clad = get_refractive_index('sio2', wavelength);
    n_box    = n_clad;  % refractive index of the substrate

    mesh_res   = round(w_sim/wavelength*100);   % approximate triangles along the longest side

    [nodes, elems, epsilon_r, regions] = build_soi_mesh( ...
        wg_width, wg_height, h_clad, h_box, w_sim, ...
        n_core, n_clad, n_box, mesh_res);
    modes = compute_modes(nodes, elems, epsilon_r, wavelength, 'num_modes', num_modes, 'mu_r', 1.0);
    tmp = zeros(1,5);

    tmp(1) = wavelength;
    tmp(2) = wg_width;
    tmp(3) = wg_height;
    tmp(4) = real(modes(1).n_eff);
    tmp(5) = real(modes(2).n_eff);

    testing_data(ii,:) = tmp;
end

save('testing_data.mat','testing_data')


figure(25); clf;
scatter3(testing_data(:,2), ...
         testing_data(:,3), ...
         testing_data(:,4), ...
         50, ...                 % marker size
         testing_data(:,4), ...  % color value
         'filled');

grid on;
xlabel('width (\mum)');
ylabel('height (\mum)');
zlabel('{\it{n}}_{eff,TE}');

colorbar;
colormap(flipud(hot));