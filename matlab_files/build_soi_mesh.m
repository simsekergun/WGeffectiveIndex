%% ============================================================
%  MESH GENERATION  (built-in, no external toolbox required)
%% ============================================================
function [nodes, elems, epsilon_r, regions] = build_soi_mesh( ...
        w_core, h_core, h_clad, h_box, w_sim, ...
        n_core, n_clad, n_box, mesh_res)
% Generate a constrained Delaunay mesh for the SOI cross-section.
%   box  : y in [-h_box, 0],     x in [-w_sim/2, w_sim/2]
%   clad : y in [0, h_clad],     x in [-w_sim/2, w_sim/2]
%   core : y in [0, h_core],     x in [-w_core/2, w_core/2]

    xL = -w_sim/2;  xR =  w_sim/2;
    yB = -h_box;    yT =  h_clad;
    xCL = -w_core/2; xCR = w_core/2;

    % resolution along each boundary (minimum 4)
    dx   = w_sim / mesh_res;
    nx   = mesh_res;
    nyb  = max(round(h_box   / dx), 4);
    nyc  = max(round(h_clad  / dx), 4);
    nxk  = max(round(w_core  / dx), 4);
    nyk  = max(round(h_core  / dx), 4);

    % ---- seed points on every interface ----
    xs   = linspace(xL,  xR,  nx)';
    xk   = linspace(xCL, xCR, nxk)';
    ysb  = linspace(yB,  0,   nyb)';
    ysc  = linspace(0,   yT,  nyc)';
    ysk  = linspace(0,   h_core, nyk)';
    ysLR = unique([ysb; ysc]);

    pts = [ xs,              repmat(yB,nx,1);      % outer bottom
            xs,              repmat(yT,nx,1);      % outer top
            repmat(xL,numel(ysLR),1), ysLR;        % outer left
            repmat(xR,numel(ysLR),1), ysLR;        % outer right
            xs,              zeros(nx,1);          % y=0 full interface
            xk,              zeros(nxk,1);         % core bottom
            xk,              repmat(h_core,nxk,1); % core top
            repmat(xCL,nyk,1), ysk;                % core left wall
            repmat(xCR,nyk,1), ysk; ];             % core right wall

    % interior fill points (avoid exact boundary overlap)
    m = dx * 0.4;
    [gx,gy] = meshgrid(linspace(xL+m, xR-m, max(round(nx*0.7),4)), ...
                        linspace(yB+m, yT-m, max(round((nyb+nyc)*0.5),4)));
    pts = [pts; gx(:), gy(:)];

    % deduplicate
    pts = uniquetol(pts, 1e-9, 'ByRows', true);

    % ---- build constraint edges along the outer boundary ----
    % Identify which points lie on each outer edge, sort them, then
    % chain consecutive pairs into constraint segments.
    snap = dx * 0.45;

    idx_bot = find(abs(pts(:,2) - yB) < snap);
    [~,o]   = sort(pts(idx_bot,1));  idx_bot = idx_bot(o);

    idx_rgt = find(abs(pts(:,1) - xR) < snap);
    [~,o]   = sort(pts(idx_rgt,2));  idx_rgt = idx_rgt(o);

    idx_top = find(abs(pts(:,2) - yT) < snap);
    [~,o]   = sort(pts(idx_top,1),'descend'); idx_top = idx_top(o);

    idx_lft = find(abs(pts(:,1) - xL) < snap);
    [~,o]   = sort(pts(idx_lft,2),'descend'); idx_lft = idx_lft(o);

    % walk around: bottom -> right -> top (reversed) -> left (reversed)
    loop = [idx_bot; idx_rgt; idx_top; idx_lft];

    % remove consecutive duplicates
    loop = loop([true; loop(2:end) ~= loop(1:end-1)]);

    % build edge pairs
    C = [loop, [loop(2:end); loop(1)]];
    C = unique(sort(C, 2), 'rows');   % canonical undirected edges

    % ---- triangulate ----
    try
        DT = delaunayTriangulation(pts, C);
    catch ME
        warning('Constrained DT failed (%s), falling back.', ME.message);
        DT = delaunayTriangulation(pts);
    end

    nodes_all = DT.Points;
    elems_all = DT.ConnectivityList;

    % ---- keep only elements whose centroid is inside the domain ----
    cx_all = (nodes_all(elems_all(:,1),1) + ...
              nodes_all(elems_all(:,2),1) + ...
              nodes_all(elems_all(:,3),1)) / 3;
    cy_all = (nodes_all(elems_all(:,1),2) + ...
              nodes_all(elems_all(:,2),2) + ...
              nodes_all(elems_all(:,3),2)) / 3;

    keep = cx_all >= xL-1e-10 & cx_all <= xR+1e-10 & ...
           cy_all >= yB-1e-10 & cy_all <= yT+1e-10;
    elems_keep = elems_all(keep, :);   % integer rows — no bounds issue

    % ---- compact node list ----
    used      = unique(elems_keep(:));
    node_map  = zeros(size(nodes_all,1), 1);
    node_map(used) = 1:numel(used);
    nodes = nodes_all(used, :);
    elems = node_map(elems_keep);

    % ---- assign permittivity ----
    Ne  = size(elems,1);
    cx  = (nodes(elems(:,1),1)+nodes(elems(:,2),1)+nodes(elems(:,3),1))/3;
    cy  = (nodes(elems(:,1),2)+nodes(elems(:,2),2)+nodes(elems(:,3),2))/3;

    epsilon_r           = ones(Ne,1) * n_clad^2;
    regions             = repmat({'clad'}, Ne, 1);

    ib = cy < -1e-10;
    epsilon_r(ib)  = n_box^2;
    regions(ib)    = {'box'};

    ic = cx >= xCL-1e-10 & cx <= xCR+1e-10 & cy >= -1e-10 & cy <= h_core+1e-10;
    epsilon_r(ic)  = n_core^2;
    regions(ic)    = {'core'};

    fprintf('  Mesh: %d nodes, %d elems  (core=%d, box=%d, clad=%d)\n', ...
        size(nodes,1), Ne, sum(ic), sum(ib), sum(~ic & ~ib));
end