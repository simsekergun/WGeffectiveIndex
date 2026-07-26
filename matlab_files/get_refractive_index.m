function n_eff = get_refractive_index(material, lambda, varargin)
% Unified refractive index calculator
%
% Inputs:
%   material : 'SiO2', 'Si3N4', 'LN', 'LNdoped'
%   lambda   : wavelength (microns preferred)
%
% Optional inputs (for LN):
%   mode     : 'even' or 'odd'
%   Ez       : electric field (V/m), default = 0
%
% Output:
%   n_eff    : effective refractive index

%% ------------------ Parse optional inputs ------------------
mode = 'even';   % default
Ez = 0;

if nargin >= 3
    mode = varargin{1};
end
if nargin >= 4
    Ez = varargin{2};
end

%% ------------------ Ensure lambda in microns ------------------
if lambda < 0.5
    lambda = lambda * 1e6;
end

%% ------------------ Material selection ------------------
switch lower(material)

    %% ========= SiO2 =========
    case 'sio2'
        n = sqrt(1 ...
            + 0.6961663./(1-(0.0684043./lambda).^2) ...
            + 0.4079426./(1-(0.1162414./lambda).^2) ...
            + 0.8974794./(1-(9.896161./lambda).^2));
        
        n_eff = n;

    %% ========= Si3N4 (Lipson) =========
    case 'si3n4'
        lambda_nm = lambda * 1e3;
        lambda2 = lambda_nm.^2;

        B1 = 3.0249;   C1 = 135.3406^2;
        B2 = 40314;    C2 = 1239842^2;

        epsr = 1 ...
            + (B1 .* lambda2) ./ (lambda2 - C1) ...
            + (B2 .* lambda2) ./ (lambda2 - C2);

        n_eff = sqrt(epsr);

    %% ========= LiNbO3 (undoped, EO included) =========
    case 'ln'
        [ne, no] = ln_sellmeier(lambda);

        % Electro-optic correction
        gamma13 = 10e-12;
        gamma33 = 33e-12;

        no = no + 0.5 * no.^3 * gamma13 .* Ez;
        ne = ne + 0.5 * ne.^3 * gamma33 .* Ez;

        n_eff = select_mode(ne, no, mode);

    %% ========= LiNbO3 doped =========
    case 'lndoped'
        [ne, no] = ln_doped_sellmeier(lambda);

        % (Typically no EO term included unless you want to add it)
        n_eff = select_mode(ne, no, mode);

    otherwise
        error('Unknown material. Choose: SiO2, Si3N4, LN, LNdoped');
end

end

%% ================= Helper: Mode selection =================
function n_eff = select_mode(ne, no, mode)

switch lower(mode)
    case 'even'
        n_eff = ne;   % extraordinary
    case 'odd'
        n_eff = no;   % ordinary
    otherwise
        error('Mode must be "even" or "odd"');
end

end

%% ================= Helper: LN Sellmeier =================
function [ne, no] = ln_sellmeier(lambda)

A = [2.9804 2.6734];
B = [0.02047 0.01764];
C = [0.5981 1.2290];
D = [0.0666 0.05914];
E = [8.9543 12.614];
F = [416.08 474.6];

% extraordinary
i = 1;
ne = sqrt(1 ...
    + A(i)*lambda.^2./(lambda.^2-B(i)) ...
    + C(i)*lambda.^2./(lambda.^2-D(i)) ...
    + E(i)*lambda.^2./(lambda.^2-F(i)));

% ordinary
i = 2;
no = sqrt(1 ...
    + A(i)*lambda.^2./(lambda.^2-B(i)) ...
    + C(i)*lambda.^2./(lambda.^2-D(i)) ...
    + E(i)*lambda.^2./(lambda.^2-F(i)));

end

%% ================= Helper: LN doped =================
function [ne, no] = ln_doped_sellmeier(lambda)

A = [2.4272 2.2454];
B = [0.01478 0.01242];
C = [1.4617 1.3005];
D = [0.05612 0.05313];
E = [9.6536 6.8972];
F = [371.216 331.33];

% extraordinary
i = 1;
ne = sqrt(1 ...
    + A(i)*lambda.^2./(lambda.^2-B(i)) ...
    + C(i)*lambda.^2./(lambda.^2-D(i)) ...
    + E(i)*lambda.^2./(lambda.^2-F(i)));

% ordinary
i = 2;
no = sqrt(1 ...
    + A(i)*lambda.^2./(lambda.^2-B(i)) ...
    + C(i)*lambda.^2./(lambda.^2-D(i)) ...
    + E(i)*lambda.^2./(lambda.^2-F(i)));

end