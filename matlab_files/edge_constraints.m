% helper: find points near a parametric line segment, sorted by t in [0,1]
function C_edge = edge_constraints(pts, p0, p1, snap_tol)
d  = p1 - p0;
L  = norm(d);
u  = d / L;
% project every point onto the line
t  = (pts - p0) .* u;           % Nx2, want row dot products
t  = t(:,1)*u(1) + t(:,2)*u(2); % scalar projection (Nx1) — note: u is 1x2
% perpendicular distance
perp = abs((pts(:,1)-p0(1))*(-u(2)) + (pts(:,2)-p0(2))*u(1));
on   = perp < snap_tol & t >= -snap_tol & t <= L+snap_tol;
idx  = find(on);
[~, o] = sort(t(idx));
idx  = idx(o);
C_edge = [idx(1:end-1), idx(2:end)];
end