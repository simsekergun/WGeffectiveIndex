%% ============================================================
%  POST-PROCESSING: TE fraction
%% ============================================================
function te_frac = compute_te_fraction(Et_dof, Ez_dof, elems, nodes, ...
                                        elem2edge, edge_sign, epsilon_r)
% Approximate TE fraction as  integral(|Et_x|^2) / integral(|Et|^2)
% evaluated by averaging over elements.

    Ne   = size(elems,1);
    local_pairs = [2 3; 3 1; 1 2];

    sum_x = 0; sum_y = 0;
    [qp, qw] = tri_quadrature();
    nq = size(qp,1);

    for ie = 1:Ne
        nds  = elems(ie,:);
        xy   = nodes(nds,:);
        edg  = elem2edge(ie,:);
        sgn  = edge_sign(ie,:);
        [~, detJ, grads] = element_geometry(xy);
        area = abs(detJ)/2;

        edge_len = zeros(1,3);
        for k = 1:3
            ni = local_pairs(k,1); nj = local_pairs(k,2);
            edge_len(k) = norm(xy(nj,:)-xy(ni,:));
        end

        for q = 1:nq
            xi = qp(q,1); eta = qp(q,2);
            lam = [1-xi-eta, xi, eta];
            w_q = qw(q) * 2 * area;

            Ex_q = 0; Ey_q = 0;
            for k = 1:3
                ni = local_pairs(k,1); nj = local_pairs(k,2);
                Wk = sgn(k)*edge_len(k)*(lam(ni)*grads(nj,:)-lam(nj)*grads(ni,:));
                Ex_q = Ex_q + Et_dof(edg(k)) * Wk(1);
                Ey_q = Ey_q + Et_dof(edg(k)) * Wk(2);
            end
            sum_x = sum_x + w_q * abs(Ex_q)^2;
            sum_y = sum_y + w_q * abs(Ey_q)^2;
        end
    end

    te_frac = real(sum_x / (sum_x + sum_y + eps));
end


