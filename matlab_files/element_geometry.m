function [B_mat, detJ, grads] = element_geometry(xy)
% Compute Jacobian and nodal gradients for a linear triangle.
% xy = [3x2] node coordinates.
    x = xy(:,1); y = xy(:,2);
    % Jacobian of affine map from reference to physical
    B_mat = [x(2)-x(1), x(3)-x(1); y(2)-y(1), y(3)-y(1)];
    detJ  = det(B_mat);
    % Gradients of barycentric coordinates (rows = nodes)
    area  = abs(detJ)/2;
    grads = zeros(3,2);
    grads(1,:) = [y(2)-y(3), x(3)-x(2)] / (2*area);
    grads(2,:) = [y(3)-y(1), x(1)-x(3)] / (2*area);
    grads(3,:) = [y(1)-y(2), x(2)-x(1)] / (2*area);
end


