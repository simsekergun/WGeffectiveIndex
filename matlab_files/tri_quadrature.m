function [qp, qw] = tri_quadrature()
% 3-point Gaussian quadrature on reference triangle (0,0)-(1,0)-(0,1)
% Exact for polynomials up to degree 2.
    qp = [1/6 1/6; 2/3 1/6; 1/6 2/3];
    qw = [1/6; 1/6; 1/6];
end

