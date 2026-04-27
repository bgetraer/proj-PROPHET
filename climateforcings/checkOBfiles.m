S=dir(fullfile(pwd,'data'));
S=S(contains({S.name},'2013') | contains({S.name},'2012'));

Nz=69;
Nx=250;
Ny=420;
Nrec=12;

for i=1:numel(S)
	fname=fullfile(S(i).folder,S(i).name);
	if contains(S(i).name,'obs')
		Nh=Nx;
	elseif contains(S(i).name,'obw')
		Nh=Ny;
	end
	A=binread(fname,8,Nh,Nz,Nrec);
	if any(isnan(A(:)))
		warning('NaN found in %s', S(i).name);
	end
end
