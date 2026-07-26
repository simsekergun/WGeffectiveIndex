clear;

lambdas= (750:25:1650)*1e-3;  %  um
wg_widths = (750:25:3000)*1e-3;  %  um
wg_heights = (600:25:1000)*1e-3;  %  um

training_data = zeros(length(lambdas)*length(wg_widths)*length(wg_heights), 5);
counter = 0;

num_modes  = 2;
compute_overlaps = 0; % set it to 1, if you want to calculate the mode overlap matrix
plot_mesh = 0;
h_clad   = 3;   % cladding height above core
h_box    = 3.00;   % buried-oxide thickness
w_sim    = 6.00;   % total simulation width


for il = 1:length(lambdas)
    wavelength = lambdas(il);
    for iw = 1:length(wg_widths)
        wg_width =wg_widths(iw);
        for ih = 1:length(wg_heights)
            wg_height = wg_heights(ih);
            counter = counter+1;
            training_data(counter,1) = wavelength;
            training_data(counter,2) = wg_width;
            training_data(counter,3) = wg_height;

            n_core = get_refractive_index('si3n4', wavelength);
            n_clad = get_refractive_index('sio2', wavelength);
            n_box    = n_clad;  % refractive index of the substrate

            mesh_res   = round(w_sim/wavelength*100);   % approximate triangles along the longest side

            [nodes, elems, epsilon_r, regions] = build_soi_mesh( ...
                wg_width, wg_height, h_clad, h_box, w_sim, ...
                n_core, n_clad, n_box, mesh_res);
            modes = compute_modes(nodes, elems, epsilon_r, wavelength, 'num_modes', num_modes, 'mu_r', 1.0);
            training_data(counter,4) = real(modes(1).n_eff);
            training_data(counter,5) = real(modes(2).n_eff);
            disp([int2str(counter) ' out of ' int2str(length(training_data))]);
        end
    end
end

save('training_data.mat','training_data')