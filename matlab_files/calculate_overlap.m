%% ============================================================
%  OVERLAP INTEGRAL  (Eq. 2 in reference)
%  overlap = 0.5 * int (E_i* x H_j + E_j x H_i*) . z  dA
%  Approximated here as  0.5 * int Et_i* . Et_j  dA
%  (proportional; full H-field reconstruction omitted for brevity)
%% ============================================================
function val = calculate_overlap(mode_i, mode_j)
% Simplified overlap using the transverse E-field inner product.
% For a proper full-vectorial overlap the H-field should be reconstructed
% from the Maxwell curl equations; this approximation suffices for
% orthogonality checks on co-propagating modes.

    elems     = mode_i.elems;
    nodes     = mode_i.nodes;
    elem2edge = mode_i.elem2edge;
    edge_sign = mode_i.edge_sign;

    Ne = size(elems,1);
    local_pairs = [2 3; 3 1; 1 2];
    [qp, qw] = tri_quadrature();
    nq = size(qp,1);

    val = 0;
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

            Ei = [0;0]; Ej = [0;0];
            for k = 1:3
                ni = local_pairs(k,1); nj = local_pairs(k,2);
                Wk = (sgn(k)*edge_len(k)*(lam(ni)*grads(nj,:)-lam(nj)*grads(ni,:)))';
                Ei = Ei + mode_i.Et_dof(edg(k)) * Wk;
                Ej = Ej + mode_j.Et_dof(edg(k)) * Wk;
            end
            val = val + w_q * (conj(Ei)' * Ej);
        end
    end
    val = 0.5 * val;
end


