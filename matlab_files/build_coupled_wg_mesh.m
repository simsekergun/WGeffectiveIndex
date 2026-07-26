function [nodes, elems, epsilon_r, regions] = build_coupled_wg_mesh( ...
    w_core, h_core, gap, ...
    h_clad, h_box, w_sim, ...
    n_core, n_clad, n_box, ...
    mesh_res)
%BUILD_COUPLED_WG_MESH  Triangular mesh for two rectangular waveguide cores
%   placed side-by-side (bus + ring straight-section) in a common cladding.
%
%   Geometry (x–y cross-section, origin at horizontal centre / BOX bottom):
%
%       |←————————————— w_sim ————————————————→|
%       ┌────────────────────────────────────────┐  ─┐
%       │              n_clad                    │   │ h_clad
%       │   ┌────────┐   gap   ┌────────┐        │  ─┤
%       │   │ n_core │◄───────►│ n_core │        │   │ h_core
%       │   └────────┘         └────────┘        │  ─┤
%       │              n_box                     │   │ h_box
%       └────────────────────────────────────────┘  ─┘
%
%   INPUTS  (all lengths in µm)
%   w_core, h_core   waveguide cross-section
%   gap              edge-to-edge gap between the two cores
%   h_clad           cladding thickness above cores
%   h_box            buried-oxide thickness below cores
%   w_sim            total simulation width  (must fit both cores + gap + margins)
%   n_core, n_clad, n_box   refractive indices
%   mesh_res         approx. elements along the longest dimension (coarse ok)
%
%   OUTPUTS  (same convention as build_soi_mesh / compute_modes)
%   nodes     [Nn×2]  node coordinates [µm]
%   elems     [Ne×3]  triangle connectivity (1-based node indices)
%   epsilon_r [Ne×1]  relative permittivity per element
%   regions   struct  geometry bookkeeping

% ── Domain extents ─────────────────────────────────────────────────────────
x_left  = -w_sim/2;
x_right =  w_sim/2;
y_bot   =  0;
y_top   =  h_box + h_core + h_clad;

y_core_bot = h_box;
y_core_top = h_box + h_core;

% Core edges (centred, gap is edge-to-edge)
x1_r =  -gap/2;          % right edge of left core  (core 1 = bus)
x1_l =  x1_r - w_core;   % left  edge of left core
x2_l =   gap/2;           % left  edge of right core (core 2 = ring)
x2_r =  x2_l + w_core;   % right edge of right core

% ── Triangle size targets ──────────────────────────────────────────────────
h_bg   = max(w_sim, y_top) / max(mesh_res, 1);  % background cell size
h_fine = min(w_core, h_core) / 8;               % inside cores
h_gap  = max(gap / 6, h_fine / 2);              % inside gap region

% ── Generate a stratified point cloud ─────────────────────────────────────
pts = point_cloud(x_left, x_right, y_bot, y_top, ...
                  x1_l, x1_r, x2_l, x2_r, ...
                  y_core_bot, y_core_top, ...
                  h_bg, h_fine, h_gap);

% ── Constrained Delaunay on the point cloud ────────────────────────────────
DT    = delaunayTriangulation(pts(:,1), pts(:,2));
tri   = DT.ConnectivityList;   % [Nt×3]
nd    = DT.Points;             % [Nn×2]

% ── Keep only triangles whose CENTROID lies inside the domain ──────────────
%    (centroid is well-defined for every valid triangle in DT.Points)
cx = (nd(tri(:,1),1) + nd(tri(:,2),1) + nd(tri(:,3),1)) / 3;
cy = (nd(tri(:,1),2) + nd(tri(:,2),2) + nd(tri(:,3),2)) / 3;

inside = cx >= x_left  & cx <= x_right & ...
         cy >= y_bot   & cy <= y_top;

tri = tri(inside, :);   % keep interior triangles only

% ── Remove unused nodes and re-index ──────────────────────────────────────
used_nodes = unique(tri(:));
new_index  = zeros(size(nd, 1), 1);
new_index(used_nodes) = 1:numel(used_nodes);

nodes = nd(used_nodes, :);            % [Nn_new × 2]
elems = new_index(tri);               % [Ne × 3]  re-indexed

% ── Assign permittivity from centroid ─────────────────────────────────────
cx = (nodes(elems(:,1),1) + nodes(elems(:,2),1) + nodes(elems(:,3),1)) / 3;
cy = (nodes(elems(:,1),2) + nodes(elems(:,2),2) + nodes(elems(:,3),2)) / 3;

in_c1 = cx >= x1_l & cx <= x1_r & cy >= y_core_bot & cy <= y_core_top;
in_c2 = cx >= x2_l & cx <= x2_r & cy >= y_core_bot & cy <= y_core_top;
in_bx = cy <  y_core_bot;

epsilon_r           = n_clad^2 * ones(size(elems,1), 1);
epsilon_r(in_bx)    = n_box^2;
epsilon_r(in_c1)    = n_core^2;
epsilon_r(in_c2)    = n_core^2;

% ── Regions struct ────────────────────────────────────────────────────────
regions.x_core1    = [x1_l, x1_r];
regions.x_core2    = [x2_l, x2_r];
regions.y_core     = [y_core_bot, y_core_top];
regions.gap        = gap;
regions.x_lim      = [x_left, x_right];
regions.y_lim      = [y_bot,  y_top];

end % main function

%% ── Helper: stratified point cloud ───────────────────────────────────────
function pts = point_cloud(xl, xr, yb, yt, ...
                            c1l, c1r, c2l, c2r, ...
                            ycb, yct, hg, hf, hgap)

% Background grid (cladding + BOX)
[Xg, Yg] = ndgrid(xl:hg:xr, yb:hg:yt);
pts = [Xg(:), Yg(:)];

% Fine grid inside core 1
[Xc, Yc] = ndgrid(c1l:hf:c1r, ycb:hf:yct);
pts = [pts; Xc(:), Yc(:)];

% Fine grid inside core 2
[Xc, Yc] = ndgrid(c2l:hf:c2r, ycb:hf:yct);
pts = [pts; Xc(:), Yc(:)];

% Very fine grid inside the gap
if c2l > c1r
    [Xg2, Yg2] = ndgrid(c1r:hgap:c2l, ycb:hgap:yct);
    pts = [pts; Xg2(:), Yg2(:)];
end

% Domain boundary corners — ensures Delaunay fills entire rectangle
pts = [pts;
       xl, yb;  xr, yb;  xr, yt;  xl, yt;
       xl, ycb; xr, ycb; xl, yct; xr, yct];

% Core 1 corners
pts = [pts;
       c1l, ycb; c1r, ycb; c1r, yct; c1l, yct];

% Core 2 corners
pts = [pts;
       c2l, ycb; c2r, ycb; c2r, yct; c2l, yct];

% Gap midline points (help triangulation resolve the gap)
if c2l > c1r
    x_mid = (c1r + c2l)/2;
    pts   = [pts; x_mid, ycb; x_mid, yct; x_mid, (ycb+yct)/2];
end

% De-duplicate (snap to half-gap grid to avoid near-coincident nodes)
snap  = hgap / 10;
pts   = round(pts / snap) * snap;
pts   = unique(pts, 'rows');

end