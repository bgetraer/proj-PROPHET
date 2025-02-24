function ptop_prime = calculateSHELFICEloadAnomaly(rhoConst)
% This function calculates the pressure load anomaly at the top of the
% water column. It is assumed that the ice shelf geometry is known, and 
% the pressure load anomaly is needed at the ice ocean interface.
%
% The pressure at the top of the ocean column is equal to:
% ptop = g sum[rho_star(k) delF(k)] + P_atm
% where:
% k is the layer index
% rho_star is the density of the displaced water at the kth layer (kg m^-3)
% delF(k) is the thickness of the kth layer
% g is the graviational acceleration
% P_atm is atmospheric pressure which is treated as negligible.
% the sum is over the fully dry cells displaced by the ice
%
% The pressure anomaly at the top of the ocean column is equal to:
% ptop_prime = g sum[(rho_star(k) - rho_const) delF(k)]
% where:
% rho_const is the reference density defined in input/data
%
% For details see documentation at:
% https://mitgcm.readthedocs.io/en/latest/phys_pkgs/shelfice.html#shelfice-description
%
% The key components of this function are
% 1) get the k indices over which to sum
%    - requires knowing where the open/closed cells are for the floating ice
% 2) calculate rho_star for each layer
%    - requires an estimate for T and S profiles for the displaced water
%    - requires a choice for the equation of state
%    - requires an iterative estimation of pressure until convergence
%      is reached
%
% p = vector of pressure profile values at cell centers
% pcnvg = root mean squared error between subsequent estimates of p
% ptol = convergence tolerance for accepting estimate of p based on pcnvg

% constants
pa2db=1E-4; % pascal to decibar

% convergence criteria
p=zc*g*rhoConst*pa2db; % initial guess for pressure in decibar (1 db = 1E-4 kg/(s^2 m))
pcnvg=rms(p); % initialize convergence criterion
ptol=1e-13; % convergence tolerance for defining the pressure

% calculate pressure, potential, etc. due to density {{{
p=zc*g*rhoConst*pa2db; % initial guess for pressure in decibar (1 db = 1E-4 kg/(s^2 m))
pcnvg=rms(p); % initialize convergence criterion
ptol=1e-13; % convergence tolerance for defining the pressure
%i=0; % benjy's loop counter, just for debugging
while pcnvg>ptol
	p0=p;                                                    % save last pressure estimate (db)
	rho=densjmd95(S_ref,T_ref,p);                            % get new estimate of density using equation of state (kg/m^3)
	drho=rho-rhoConst;                                       % density anomaly (kg/m^3)
	phiC=cumsum(dz*g*drho/rhoConst)-(dz/2)*g*drho/rhoConst;  % cumulative gravitational potential anomaly at the cell centers (m^2/s^2)
	phiF=[0 cumsum(dz*g*drho/rhoConst)];                     % cumulative gravitational potential anomaly at ALL cell edges (m^2/s^2)
	p=rhoConst*(zc*g+phiC)*pa2db;                            % new pressure estimate (db)
	pcnvg=rms(p-p0);                                         % update convergence criterion
	%i=i+1; disp(num2str(i)); % loop counter for debugging
end
massC=rhoConst*(phiC/g+ zc); % vertically integrated mass density at cell centers (kg/m^2)
massF=rhoConst*(phiF/g+ zp); % vertically integrated mass density at ALL cell edges (kg/m^2)
