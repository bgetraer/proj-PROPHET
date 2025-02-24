zc=mit.mesh.zc;
rhoConst = 1030;
g = 9.81;
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
   rho=densjmd95(S,T,p);                            % get new estimate of density using equation of state (kg/m^3)
   drho=rho-rhoConst;                                       % density anomaly (kg/m^3)
   phiC=cumsum(dz*g*drho/rhoConst)-(dz/2)*g*drho/rhoConst;  % cumulative gravitational potential anomaly at the cell centers (m^2/s^2)
   phiF=[0 cumsum(dz*g*drho/rhoConst)];                     % cumulative gravitational potential anomaly at ALL cell edges (m^2/s^2)
   p=rhoConst*(zc*g+phiC)*pa2db;                            % new pressure estimate (db)
   pcnvg=rms(p-p0);                                         % update convergence criterion
   %i=i+1; disp(num2str(i)); % loop counter for debugging
end
