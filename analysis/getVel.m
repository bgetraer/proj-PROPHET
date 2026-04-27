
fname = 'experiments/Paris2C/RUN02/results/pickup.save.2021760000.data';
addpath('/totten_1/bgetraer/issmjpl/proj-getraer/issmxmitgcm/issmxmitgcm');
mit = loadmodel('experiments/MITgcm_initialization/Models/PROPHET_mitgcm_init_MeshInit.mat');
P = readpickup(fname,mit.mesh.Nx,mit.mesh.Ny,mit.mesh.Nz);

U_q = interp3(mit.mesh.XC - mit.mesh.delxF, mit.mesh.YC, mit.mesh.ZC, P.Uvel, mit.mesh.XC, mit.mesh.YC, mit.mesh.ZC);
V_q = interp3(mit.mesh.XC, mit.mesh.YC - mit.mesh.delyF, mit.mesh.ZC, P.Vvel, mit.mesh.XC, mit.mesh.YC, mit.mesh.ZC);

Vel_q = (U_q.^2 + V_q.^2);

