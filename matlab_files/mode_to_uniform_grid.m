function [Ex, Ey, Ez, Hx, Hy, Hz] = mode_to_uniform_grid(mode, grid_x, grid_y)
% MODE_TO_UNIFORM_GRID  Interpolate FEM mode fields onto a uniform Cartesian grid.
%
%   [Ex, Ey, Ez, Hx, Hy, Hz] = mode_to_uniform_grid(mode, grid_x, grid_y)
%
%   Inputs:
%     mode    - struct from compute_modes (single mode, not array)
%     grid_x  - 1×Nx vector of x-coordinates  [µm]
%     grid_y  - 1×Ny vector of y-coordinates  [µm]
%
%   Outputs:
%     Ex, Ey, Ez  - Ny×Nx complex electric field matrices
%     Hx, Hy, Hz  - Ny×Nx complex magnetic field matrices
%                   (NaN outside the FEM domain)

nodes    = mode.nodes;
elems    = mode.elems;
elem2edge = mode.elem2edge;
edge_sign = mode.edge_sign;

Nx = numel(grid_x);
Ny = numel(grid_y);

[GX, GY] = meshgrid(grid_x, grid_y);          % Ny×Nx
pts = [GX(:), GY(:)];                          % (Nx*Ny) × 2

Ne = size(elems, 1);
local_pairs = [2 3; 3 1; 1 2];

% ----------------------------------------------------------------
% 1.  Build element bounding boxes for fast point-in-triangle search
% ----------------------------------------------------------------
xn = nodes(:,1);  yn = nodes(:,2);

x1 = xn(elems(:,1)); x2 = xn(elems(:,2)); x3 = xn(elems(:,3));
y1 = yn(elems(:,1)); y2 = yn(elems(:,2)); y3 = yn(elems(:,3));

bbxmin = min([x1 x2 x3],[],2);
bbxmax = max([x1 x2 x3],[],2);
bbymin = min([y1 y2 y3],[],2);
bbymax = max([y1 y2 y3],[],2);

% ----------------------------------------------------------------
% 2.  For every grid point find its host triangle & barycentric coords
% ----------------------------------------------------------------
Np        = size(pts,1);
host_elem = zeros(Np,1,'int32');   % 0 = not found (outside domain)
bary      = zeros(Np,3);

for ip = 1:Np
    px = pts(ip,1);  py = pts(ip,2);

    % Candidate elements via bounding box
    cands = find(px >= bbxmin & px <= bbxmax & ...
                 py >= bbymin & py <= bbymax);

    for ic = 1:numel(cands)
        ie = cands(ic);
        ax = x1(ie); ay = y1(ie);
        bx = x2(ie); by = y2(ie);
        cx = x3(ie); cy = y3(ie);

        % Barycentric coordinates
        denom = (by-cy)*(ax-cx) + (cx-bx)*(ay-cy);
        l1 = ((by-cy)*(px-cx) + (cx-bx)*(py-cy)) / denom;
        l2 = ((cy-ay)*(px-cx) + (ax-cx)*(py-cy)) / denom;
        l3 = 1 - l1 - l2;

        tol = -1e-10;
        if l1 >= tol && l2 >= tol && l3 >= tol
            host_elem(ip) = ie;
            bary(ip,:)    = [l1, l2, l3];
            break
        end
    end
end

% ----------------------------------------------------------------
% 3.  Evaluate E-fields at every located grid point
% ----------------------------------------------------------------
Ex_v = nan(Np,1,'like',1j+0);
Ey_v = nan(Np,1,'like',1j+0);
Ez_v = nan(Np,1,'like',1j+0);

for ip = 1:Np
    ie = host_elem(ip);
    if ie == 0, continue, end

    nds = elems(ie,:);
    xy  = nodes(nds,:);
    lam = bary(ip,:);

    [~, ~, grads] = element_geometry(xy);

    edg = elem2edge(ie,:);
    sgn = edge_sign(ie,:);

    edge_len = zeros(1,3);
    for k = 1:3
        ni = local_pairs(k,1); nj = local_pairs(k,2);
        edge_len(k) = norm(xy(nj,:)-xy(ni,:));
    end

    ex = 0; ey = 0;
    for k = 1:3
        ni = local_pairs(k,1); nj = local_pairs(k,2);
        Wk = sgn(k)*edge_len(k)*(lam(ni)*grads(nj,:) - lam(nj)*grads(ni,:));
        ex = ex + mode.Et_dof(edg(k)) * Wk(1);
        ey = ey + mode.Et_dof(edg(k)) * Wk(2);
    end

    % Ez: nodal (linear Lagrange) interpolation
    ez = lam(1)*mode.Ez_dof(nds(1)) + ...
         lam(2)*mode.Ez_dof(nds(2)) + ...
         lam(3)*mode.Ez_dof(nds(3));

    Ex_v(ip) = ex;
    Ey_v(ip) = ey;
    Ez_v(ip) = ez;
end

% ----------------------------------------------------------------
% 4.  Derive H-fields from curl E = -iωμ₀ H
%     H_t from curl_z(E_t) and ∂Ez/∂x,y
%     Hz  from curl_t(E_t)
%     Uses finite differences on the uniform grid for derivatives.
% ----------------------------------------------------------------
dx = grid_x(2)-grid_x(1);
dy = grid_y(2)-grid_y(1);
k0 = mode.k0;
beta = mode.k;                % propagation constant (k0 * n_eff)
mu0  = 4*pi*1e-7;
eps0 = 8.854187817e-12;
c0   = 1/sqrt(mu0*eps0);
omega = mode.omega;
Z0    = mu0*omega;            % ωμ₀  (SI factor cancels with geometry µm→m below)

% Reshape to 2D grids
Ex2 = reshape(Ex_v, Ny, Nx);
Ey2 = reshape(Ey_v, Ny, Nx);
Ez2 = reshape(Ez_v, Ny, Nx);

% ∂/∂x, ∂/∂y with central differences (NaN-aware via simple padding)
dEz_dx = central_diff(Ez2, dx, 2);   % d/dx  → along columns
dEz_dy = central_diff(Ez2, dy, 1);   % d/dy  → along rows
dEx_dy = central_diff(Ex2, dy, 1);
dEy_dx = central_diff(Ey2, dx, 2);

% Maxwell: curl E = -iωμ₀ H  (with e^{iβz} dependence, ∂/∂z → iβ)
%
%   Hx = (1/(-iωμ₀)) * ( iβ·Ey - ∂Ez/∂y )
%   Hy = (1/(-iωμ₀)) * ( ∂Ez/∂x - iβ·Ex )
%   Hz = (1/(-iωμ₀)) * ( ∂Ey/∂x - ∂Ex/∂y )

fac = 1 / (-1i * Z0);

Hx2 = fac * (1i*beta*Ey2 - dEz_dy);
Hy2 = fac * (dEz_dx - 1i*beta*Ex2);
Hz2 = fac * (dEy_dx - dEx_dy);

% ----------------------------------------------------------------
% 5.  Output
% ----------------------------------------------------------------
Ex = Ex2;  Ey = Ey2;  Ez = Ez2;
Hx = Hx2;  Hy = Hy2;  Hz = Hz2;
end


% ----------------------------------------------------------------
%  Helper: central finite difference along dimension dim
% ----------------------------------------------------------------
function df = central_diff(f, h, dim)
df = zeros(size(f),'like',f);
if dim == 1   % d/dy
    df(2:end-1,:) = (f(3:end,:) - f(1:end-2,:)) / (2*h);
    df(1,:)       = (f(2,:)     - f(1,:))        / h;
    df(end,:)     = (f(end,:)   - f(end-1,:))    / h;
else           % d/dx
    df(:,2:end-1) = (f(:,3:end) - f(:,1:end-2)) / (2*h);
    df(:,1)       = (f(:,2)     - f(:,1))        / h;
    df(:,end)     = (f(:,end)   - f(:,end-1))    / h;
end
end