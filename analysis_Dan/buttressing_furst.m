function Kn_vertices=buttressing_furst(md,vx,vy,thickness)
%Function that plots buttressing ratio over an ice shelf, approximating 
%Gudmundsson's definition thetaN = N/N0 in the direction of flow as opposed to the normal to the grounding line.
% N is the sigma_nn
% N0 would be the stress exerted by the ocean
% theta = 1 means no buttressing
% theta close to 0 is a lot of buttressing

%Find Elements that are floating
index = md.mesh.elements;
pos_floating = find(min(md.mask.ocean_levelset(index),[],2)<0);

%Recover stresses
md=mechanicalproperties(md, vx, vy);

%Hilmar's horrible cmap
cmap = jet(80);
cmap(60-5:60+5,:) = 0;

%Allocate thetaN
thetaN  = zeros(md.mesh.numberofelements,1);

for el=1:md.mesh.numberofelements
	% normal vector is taken to be along flow
	vx_el = mean(vx(index(el,:)));
	vy_el = mean(vy(index(el,:)));
	vel_el = sqrt(vx_el^2+vy_el^2);
	n=[vx_el;vy_el]/vel_el;

	%Build sigma_nn
	tau_xx = md.results.deviatoricstress.xx(el);
	tau_yy = md.results.deviatoricstress.yy(el);
	tau_xy = md.results.deviatoricstress.xy(el);
	R = [2*tau_xx+tau_yy   tau_xy;tau_xy   2*tau_yy+tau_xx];
	N = n'*R*n;

	%Water stress only
	H = mean(thickness(index(el,:)));
	g = md.constants.g;
	rho_i = md.materials.rho_ice;
	rho_w = md.materials.rho_water;
	N0 = 0.5*g*rho_i.*(1-rho_i./rho_w).*H;

	%Plot thetaN
	thetaN(el)  = N/N0;
end

Kn = 1-thetaN;
Kn(Kn<0) = 0;
Kn(Kn>1) = 1;

Kn_vertices = averaging(md,Kn,0);
Kn_vertices(md.mask.ocean_levelset>0) = NaN;
%Kn_vertices(isnan(Kn_vertices)) = 0;
