mit=loadmodel('../experiments/MITgcm_initialization/Models/PROPHET_mitgcm_init_RuntimeOptionsCoupled.mat');
Nz = 69;
T=binread('../experiments/MITgcm_initialization/input/tref.bin',8,[1,Nz]);
S=binread('../experiments/MITgcm_initialization/input/sref.bin',8,[1,Nz]);

g = 9.81; % Gravitational acceleration (m/s^2)
rho_const = 1027; % Reference density (kg/m^3) - Example value
ptol = 1e-3; % Convergence tolerance (Pa)
delF = mit.mesh.delzF;
%[p, rho_star, pcnvg] = calc_pressure_profile(mit.mesh.delzF, T, S, rho_const, g, ptol)

%function [p, rho_star, pcnvg] = calc_pressure_profile(delF, T, S, rho_const, g, ptol)
% calculate_pressure_profile - Iteratively calculates pressure profile.
%
% This subfunction iteratively calculates the pressure profile and
% density profile until convergence is reached.
%
% Inputs:
%   delF      - Vector of layer thicknesses (m).
%   T         - Vector of temperature profile (degC).
%   S         - Vector of salinity profile (PSU).
%   rho_const - Reference density (kg/m^3).
%   g         - Gravitational acceleration (m/s^2).
%   ptol      - Convergence tolerance (Pa).
%
% Outputs:
%   p         - Row vector of pressure profile values (Pa).
%   rho_star  - Row vector of density profile values (kg/m^3).
%   pcnvg     - Root mean squared error between successive pressure estimates (Pa).
pa2db=1E-4; % pascal to decibar

maxIterations = 100; % Maximum iterations for pressure convergence
numLayers = length(T);
p = g*rho_const*mit.mesh.zc*pa2db; 
p_old = zeros(numLayers,1);
pcnvg = sqrt(mean((p - p_old).^2));

for iter = 1:maxIterations

	rho=densjmd95(S,T,p); % get new estimate of density using equation of state (kg/m^3)

    % Calculate pressure profile
    p = cumsum(g * rho .* delF);

    % Check for convergence
    pcnvg = sqrt(mean((p - p_old).^2));
    if pcnvg < ptol
        return;
    end
    p_old = p;
end
%end
