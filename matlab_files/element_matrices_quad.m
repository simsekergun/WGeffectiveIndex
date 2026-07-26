%% ============================================================
%  ELEMENT MATRICES (quadrature-based)
%% ============================================================
function [Ae, Be] = element_matrices_quad( ...
        xy, grads, area, edge_len, sgn, qp, qw, eps, mu_r, k0)
% Assemble 6x6 element matrices using numerical quadrature.
% Local DOFs: [edge1 edge2 edge3 node1 node2 node3]

    nq  = size(qp,1);
    Ae  = zeros(6,6,'like',eps+1j);
    Be  = zeros(6,6,'like',eps+1j);

    local_edges = [2 3; 3 1; 1 2];  % local node pairs

    % Nodal gradients (constant over element for P1)
    % grads = [3x2]: grad(phi_k) for k=1,2,3
    mu_t = mu_r;
    mu_z = mu_r;

    for q = 1:nq
        xi  = qp(q,1);
        eta = qp(q,2);
        w   = qw(q) * 2 * area;   % physical weight

        % Barycentric coordinates at quad point
        lam = [1-xi-eta, xi, eta];

        % Nedelec shape functions (vector, physical coords) [3x2]
        % w_k = s_k * l_k * (lam_i * grad_j - lam_j * grad_i)
        Wt = zeros(3,2);
        for k = 1:3
            ni = local_edges(k,1);
            nj = local_edges(k,2);
            Wt(k,:) = sgn(k) * edge_len(k) * ...
                      (lam(ni)*grads(nj,:) - lam(nj)*grads(ni,:));
        end

        % curl of Nedelec functions (scalar in 2D) = curl_z w_k
        % curl_z w_k = s_k * l_k * 2 * (grad_i x grad_j)
        curlW = zeros(3,1);
        for k = 1:3
            ni = local_edges(k,1);
            nj = local_edges(k,2);
            gi = grads(ni,:);
            gj = grads(nj,:);
            curlW(k) = sgn(k) * edge_len(k) * 2 * (gi(1)*gj(2) - gi(2)*gj(1));
        end

        % P1 nodal shape functions and gradients
        Phi  = lam;          % [1x3]
        GPhi = grads;        % [3x2]: GPhi(k,:) = grad phi_k

        % ---- Assemble contributions ----
        for ii = 1:6
            for jj = 1:6
                % Determine if i,j are edge (1-3) or nodal (4-6) DOFs
                if ii <= 3 && jj <= 3
                    % (A) curl-curl term:  (1/mu_z) * curl(w_i) * curl(w_j) / k0^2
                    a_cc = (1/mu_z) * curlW(ii) * curlW(jj) / k0^2;
                    % (B) mass term (transverse):  -eps_t * w_i . w_j
                    a_mm = -eps * dot(Wt(ii,:), Wt(jj,:));
                    Ae(ii,jj) = Ae(ii,jj) + w * (a_cc + a_mm);
                    % B-form: -(1/mu_t) * w_i . w_j / k0^2
                    Be(ii,jj) = Be(ii,jj) - w * (1/mu_t) * dot(Wt(ii,:), Wt(jj,:)) / k0^2;

                elseif ii <= 3 && jj > 3  % edge test, nodal trial
                    jn = jj - 3;
                    % (1/mu_t) * grad(phi_j) . w_i
                    a_gz = (1/mu_t) * dot(GPhi(jn,:), Wt(ii,:));
                    Ae(ii,jj) = Ae(ii,jj) + w * a_gz;

                elseif ii > 3 && jj <= 3  % nodal test, edge trial
                    in = ii - 3;
                    % eps_t * w_j . grad(phi_i)
                    a_et = eps * dot(Wt(jj,:), GPhi(in,:));
                    Ae(ii,jj) = Ae(ii,jj) + w * a_et;

                else  % nodal-nodal  (ii>3, jj>3)
                    in = ii - 3;
                    jn = jj - 3;
                    % -eps_z * phi_i * phi_j * k0^2
                    a_zz = -eps * Phi(in) * Phi(jn) * k0^2;
                    Ae(ii,jj) = Ae(ii,jj) + w * a_zz;
                    % B: no nodal-nodal contribution
                end
            end
        end
    end
end
