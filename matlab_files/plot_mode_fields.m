%% ============================================================
%  VISUALIZATION
%% ============================================================
function plot_mode_fields(mode, nodes, elems, ttl)
    Ne = size(elems,1);
    local_pairs = [2 3; 3 1; 1 2];
    elem2edge = mode.elem2edge;
    edge_sign = mode.edge_sign;

    % Evaluate fields at element centroids
    Ex_c = zeros(Ne,1,'like',1j);
    Ey_c = zeros(Ne,1,'like',1j);
    Ez_c = zeros(Ne,1,'like',1j);

    for ie = 1:Ne
        nds  = elems(ie,:);
        xy   = nodes(nds,:);
        edg  = elem2edge(ie,:);
        sgn  = edge_sign(ie,:);
        [~, ~, grads] = element_geometry(xy);
        area_unused = 1;

        % Centroid = (1/3, 1/3) in barycentric
        lam = [1/3, 1/3, 1/3];
        edge_len = zeros(1,3);
        for k = 1:3
            ni = local_pairs(k,1); nj = local_pairs(k,2);
            edge_len(k) = norm(xy(nj,:)-xy(ni,:));
        end

        for k = 1:3
            ni = local_pairs(k,1); nj = local_pairs(k,2);
            Wk = sgn(k)*edge_len(k)*(lam(ni)*grads(nj,:)-lam(nj)*grads(ni,:));
            Ex_c(ie) = Ex_c(ie) + mode.Et_dof(edg(k)) * Wk(1);
            Ey_c(ie) = Ey_c(ie) + mode.Et_dof(edg(k)) * Wk(2);
        end
        % Ez at centroid (average nodal values)
        Ez_c(ie) = mean(mode.Ez_dof(nds));
    end

    % Normalise
    Emax = max([max(abs(Ex_c)), max(abs(Ey_c)), max(abs(Ez_c))]);
    Ex_c = Ex_c / Emax; Ey_c = Ey_c / Emax; Ez_c = Ez_c / Emax;

    % Centroid coordinates for patch plot
    cx = mean(nodes(elems,1), 2);
    cy = mean(nodes(elems,2), 2);

    figure('Name', ttl, 'Position', [100 100 1200 380]);
    cmps = {abs(Ex_c), abs(Ey_c), abs(Ez_c)};
    lbls = {'|E_x|', '|E_y|', '|E_z|'};

    for k = 1:3
        ax = subplot(1,3,k);
        tripcolor_custom(ax, nodes, elems, cmps{k});
        title(lbls{k}, 'Interpreter','tex');
        xlabel('x  [\mum]'); ylabel('y  [\mum]');
        colorbar; axis equal tight;
    end
    sgtitle(sprintf('%s   n_{eff} = %.6f + %.2ei', ...
        ttl, real(mode.n_eff), imag(mode.n_eff)), 'Interpreter','tex');
    % print
    filename = [strrep(ttl, ' ', '_'), '.png'];
    print(gcf, filename, '-dpng', '-r300');
end


