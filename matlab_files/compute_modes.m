%% ============================================================
%  FUNCTION: compute_modes
%% ============================================================
function modes = compute_modes(nodes, elems, epsilon_r, wavelength, varargin)
% COMPUTE_MODES  Full-vectorial FEM mode solver for dielectric waveguides.
%
%   modes = compute_modes(nodes, elems, epsilon_r, wavelength)
%   modes = compute_modes(... , 'num_modes', 4, 'mu_r', 1, 'n_guess', [])
%
%   nodes      – [Nn x 2]  node coordinates (x,y)
%   elems      – [Ne x 3]  triangle connectivity (1-based)
%   epsilon_r  – [Ne x 1]  relative permittivity per element (complex OK)
%   wavelength – scalar    free-space wavelength (same units as nodes)
%
%   Returns a struct array with fields:
%     k, n_eff, omega, k0, wavelength, te_fraction, tm_fraction,
%     Et_dof, Ez_dof,            % raw DOF vectors
%     nodes, elems, epsilon_r    % mesh (for plotting / post-processing)

    %% Parse optional arguments
    p = inputParser;
    addParameter(p, 'num_modes',          4,    @isnumeric);
    addParameter(p, 'mu_r',               1.0,  @isnumeric);
    addParameter(p, 'n_guess',            [],   @isnumeric);
    addParameter(p, 'metallic_boundaries',false, @islogical);
    parse(p, varargin{:});
    opts = p.Results;

    c0     = 2.99792458e14;   % speed of light  µm/s  (matches µm mesh)
    omega  = 2*pi*c0 / wavelength;
    k0     = omega / c0;   % = 2*pi/lambda

    Nn = size(nodes,1);
    Ne = size(elems,1);

    %% Build local element data (Nedelec + Lagrange on triangles)
    % We use the lowest-order mixed element:
    %   - 3 edge (Nedelec) DOFs per triangle for E_t
    %   - 3 nodal (Lagrange P1) DOFs per triangle for E_z
    %
    % Global DOF layout:  [Et_edges (Nedge) | Ez_nodes (Nn)]
    %                      1..Nedge          Nedge+1..Nedge+Nn

    % Build edge table
    [edges, elem2edge, edge_sign] = build_edge_table(elems);
    Nedge = size(edges,1);
    Ndof  = Nedge + Nn;

    fprintf('  DOFs: %d edge (Et) + %d nodal (Ez) = %d total\n', ...
        Nedge, Nn, Ndof);

    %% Assemble global matrices A and B
    % Eigenvalue problem:  A * x = lambda * B * x
    % where lambda = beta^2 / k0^2  (= n_eff^2)
    %
    % Weak form (from Eq. (7)-(8) of the reference):
    %
    %   A(u,v) =  (1/mu_z) * curl(e_t) * curl(v_t) / k0^2
    %           - epsilon_t * e_t . v_t
    %           + (1/mu_t) * grad(e_z) . v_t
    %           + epsilon_t * e_t . grad(v_z)
    %           - epsilon_z * e_z * v_z * k0^2
    %
    %   B(u,v) = -(1/mu_t) * e_t . v_t / k0^2
    %
    % The sign convention gives eigenvalues lam = beta^2/k0^2 > 0 for
    % guided modes, and n_eff = sqrt(lam).

    % Pre-allocate sparse triplet storage (upper bound: 36 entries/elem)
    nnz_est = 36 * Ne;
    Ai = zeros(nnz_est,1,'int32');
    Aj = zeros(nnz_est,1,'int32');
    Av = zeros(nnz_est,1,'like', epsilon_r(1)+1j);
    Bi = zeros(nnz_est,1,'int32');
    Bj = zeros(nnz_est,1,'int32');
    Bv = zeros(nnz_est,1,'like', epsilon_r(1)+1j);
    cnt_A = 0;
    cnt_B = 0;

    mu_r = opts.mu_r;

    for ie = 1:Ne
        nds  = elems(ie,:);           % global node indices  [1x3]
        xy   = nodes(nds,:);          % [3x2] node coords
        eps  = epsilon_r(ie);         % scalar permittivity
        edg  = elem2edge(ie,:);       % global edge indices  [1x3]
        sgn  = edge_sign(ie,:);       % ±1 orientation signs [1x3]

        % Compute element geometry
        [B_mat, detJ, grads] = element_geometry(xy);
        area = abs(detJ) / 2;

        % ----------------------------------------------------------------
        %  Local Nedelec shape functions (lowest order, Whitney 1-forms)
        %  On the reference element the basis functions are:
        %    w_k = s_k * l_k * (lam_i * grad(lam_j) - lam_j * grad(lam_i))
        %  where edge k connects local nodes i->j, l_k = edge length,
        %  s_k = orientation sign.
        %
        %  In physical coordinates we evaluate curl and value at
        %  Gaussian quadrature points.
        % ----------------------------------------------------------------

        % Edge lengths
        edge_len = zeros(1,3);
        local_edges = [2 3; 3 1; 1 2];   % local node pairs for edges 1,2,3
        for k = 1:3
            ni = local_edges(k,1); nj = local_edges(k,2);
            edge_len(k) = norm(xy(nj,:) - xy(ni,:));
        end

        % Gaussian quadrature on triangle (3-point rule, exact for deg 2)
        [qp, qw] = tri_quadrature();   % [3x2] ref coords, [3x1] weights
        nq = size(qp,1);

        % Accumulate element matrices via quadrature
        [Ae_local, Be_local] = element_matrices_quad( ...
            xy, grads, area, edge_len, sgn, qp, qw, eps, mu_r, k0);

        % ----------------------------------------------------------------
        %  Scatter local (6x6) matrices into global triplets
        %  Local DOF ordering: [edge1 edge2 edge3 | node1 node2 node3]
        % ----------------------------------------------------------------
        local_dofs = [edg, (nds + Nedge)];   % 1x6  global DOF indices

        for ii = 1:6
            for jj = 1:6
                gi = local_dofs(ii);
                gj = local_dofs(jj);
                aval = Ae_local(ii,jj);
                bval = Be_local(ii,jj);

                if aval ~= 0
                    cnt_A = cnt_A + 1;
                    Ai(cnt_A) = gi;
                    Aj(cnt_A) = gj;
                    Av(cnt_A) = aval;
                end
                if bval ~= 0
                    cnt_B = cnt_B + 1;
                    Bi(cnt_B) = gi;
                    Bj(cnt_B) = gj;
                    Bv(cnt_B) = bval;
                end
            end
        end
    end

    A = sparse(double(Ai(1:cnt_A)), double(Aj(1:cnt_A)), Av(1:cnt_A), Ndof, Ndof);
    B = sparse(double(Bi(1:cnt_B)), double(Bj(1:cnt_B)), Bv(1:cnt_B), Ndof, Ndof);

    %% Apply boundary conditions (PEC / Dirichlet on boundary edges)
    if opts.metallic_boundaries
        disp('line 229')
        bnd_dofs = boundary_edge_dofs(edges, elems, Nn);
    else
        % PML-free open boundary: just fix the Ez DOFs on outer boundary
        bnd_nodes = boundary_nodes(elems, Nn);
        bnd_dofs  = bnd_nodes + Nedge;
    end

    free_dofs = setdiff(1:Ndof, bnd_dofs);
    Af = A(free_dofs, free_dofs);
    Bf = B(free_dofs, free_dofs);

    %% Shift-invert sparse eigensolver
    if isempty(opts.n_guess)
        sigma = k0^2 * max(real(epsilon_r)) * 1.1;
    else
        sigma = k0^2 * opts.n_guess^2;
    end

    fprintf('  Running eigs (k=%d, sigma=%.4f) ...\n', opts.num_modes, sigma);
    opts_eigs.tol     = 1e-10;
    opts_eigs.maxit   = 500;
    opts_eigs.issym   = false;
    opts_eigs.isreal  = false;

    try
        [Vf, D] = eigs(-Af, -Bf, opts.num_modes, sigma, opts_eigs);
    catch ME
        warning('eigs failed: %s\nRetrying without shift.', ME.message);
        [Vf, D] = eigs(-Af, -Bf, opts.num_modes, 'lm', opts_eigs);
    end

    lams = diag(D);    % eigenvalues = n_eff^2

    %% Sort by real part of n_eff (descending = most confined first)
    [~, ord] = sort(real(sqrt(lams)), 'descend');
    lams = lams(ord);
    Vf   = Vf(:,ord);

    %% Reconstruct full DOF vectors & build mode structs
    modes = struct();
    for m = 1:opts.num_modes
        x_full = zeros(Ndof,1,'like',1j);
        x_full(free_dofs) = Vf(:,m);

        Et_dof = x_full(1:Nedge);
        Ez_dof = x_full(Nedge+1:end);

        lam   = lams(m);
        % beta  = sqrt(lam) * k0;   % propagation constant
        % n_eff = sqrt(lam);

        beta  = sqrt(lam);  % propagation constant (1/μm)
        n_eff = beta / k0;  % effective index

        % TE/TM fraction via energy in Et vs Ez
        [te_frac] = compute_te_fraction(Et_dof, Ez_dof, elems, ...
                                         nodes, elem2edge, edge_sign, epsilon_r);

        modes(m).k          = beta;
        modes(m).n_eff      = n_eff;
        modes(m).lambda     = wavelength;
        modes(m).omega      = omega;
        modes(m).k0         = k0;
        modes(m).te_fraction = te_frac;
        modes(m).tm_fraction = 1 - te_frac;
        modes(m).Et_dof     = Et_dof;
        modes(m).Ez_dof     = Ez_dof;
        modes(m).edges      = edges;
        modes(m).elem2edge  = elem2edge;
        modes(m).edge_sign  = edge_sign;
        modes(m).nodes      = nodes;
        modes(m).elems      = elems;
        modes(m).epsilon_r  = epsilon_r;
        modes(m).Nedge      = Nedge;
    end
    format long
    n_eff_list = [modes.n_eff].'
end


