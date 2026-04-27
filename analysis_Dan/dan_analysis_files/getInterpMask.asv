% LOAD SHAPE FILES AND INTERPOLATION MESH FROM DAN
filename = 'IceBoundaries_Y2014-2016_Antarctica/IceBoundaries_Y2014-2016_Antarctica.shp';
A = shaperead(filename);
is_grounded = strcmp({A.TYPE},"GR"); % index of grounded ice domains
is_floating = strcmp({A.TYPE},"FL"); % index of floating ice domains

filename = 'for_averaging.mat';
B = load(filename);
[B.X, B.Y] = meshgrid(B.x_mesh,B.y_mesh); % mesh grid of the x and y cell boundaries
B.XC = B.X(1:end-1,1:end-1) + diff(B.X(1:end-1,:),[],2)/2; % mesh grid of the x cell centers
B.YC = B.Y(1:end-1,1:end-1) + diff(B.Y(:,1:end-1),[],1)/2; % mesh grid of the y cell centers
B.mask = B.mask'==1;

% PIG — combination of grounded Pine_Island and floating Pine_Island
name = "Pine_Island";
PIG_GR = A(strcmp({A.NAME},name) & is_grounded);
PIG_FL = A(strcmp({A.NAME},name) & is_floating);

PIG_polyvec_gr = [polyshape(PIG_GR.X,PIG_GR.Y)]; % region grounded polyvec 
PIG_polyvec_fl = [polyshape(PIG_FL.X,PIG_FL.Y)]; % region floating polyvec 
[PIG_mask, PIG_shape] = mask_region(PIG_polyvec_gr,PIG_polyvec_fl,B.XC,B.YC,B.mask); % extract region shape and mask

% THW — combination of grounded Thwaites, floating Thwaites, and grounded Haynes
name = "Thwaites";
THWAITES_GR = A(strcmp({A.NAME},name) & is_grounded);
THWAITES_FL = A(strcmp({A.NAME},name) & is_floating);
name = "Haynes";
HAYNES_GR = A(strcmp({A.NAME},name) & is_grounded);

THW_polyvec_gr = [polyshape(THWAITES_GR.X,THWAITES_GR.Y),... 
    polyshape(HAYNES_GR.X,HAYNES_GR.Y)]; % region grounded polyvec 
THW_polyvec_fl = [polyshape(THWAITES_FL.X,THWAITES_FL.Y)]; % region floating polyvec 

[THW_mask, THW_shape] = mask_region(THW_polyvec_gr,THW_polyvec_fl,B.XC,B.YC,B.mask); % extract region shape and mask

% SMITH - combination of grounded Kohler, grounded Pope, grounded Smith, floating Dotson, and floating Crosson.
name = "Kohler";
KOHLER_GR = A(strcmp({A.NAME},name) & is_grounded);
name = "Pope";
POPE_GR = A(strcmp({A.NAME},name) & is_grounded);
name = "Smith";
SMITH_GR = A(strcmp({A.NAME},name) & is_grounded);
name = "Dotson";
DOTSON_FL = A(strcmp({A.NAME},name) & is_floating);
name = "Crosson";
CROSSON_FL = A(strcmp({A.NAME},name) & is_floating);

SMITH_polyvec_gr = [polyshape(KOHLER_GR.X,KOHLER_GR.Y),...
    polyshape(POPE_GR.X,POPE_GR.Y),...
    polyshape(SMITH_GR.X,SMITH_GR.Y)]; % region grounded polyvec 
SMITH_polyvec_fl = [polyshape(DOTSON_FL.X,DOTSON_FL.Y),...
    polyshape(CROSSON_FL.X,CROSSON_FL.Y)];  % region floating polyvec 
[SMITH_mask, SMITH_shape] = mask_region(SMITH_polyvec_gr,SMITH_polyvec_fl,B.XC,B.YC,B.mask); % extract region shape and mask

% PLOT
figure(1);clf;hold on;
% plot(PIG_shape,'FaceColor','b','EdgeColor','none')
plot(B.XC(PIG_mask),B.YC(PIG_mask),'.b')
contour(B.XC,B.YC,PIG_mask,[0.5,0.5],'EdgeColor','k','LineWidth',2)
plot(PIG_shape,'FaceColor','none','EdgeColor','k','LineWidth',2)
% plot(THW_shape,'FaceColor','y','EdgeColor','none')
plot(B.XC(THW_mask),B.YC(THW_mask),'.g')
contour(B.XC,B.YC,THW_mask,[0.5,0.5],'EdgeColor','k','LineWidth',2)
plot(THW_shape,'FaceColor','none','EdgeColor','k','LineWidth',2)
% plot(SMITH_shape,'FaceColor','r','EdgeColor','none')
plot(B.XC(SMITH_mask),B.YC(SMITH_mask),'.r')
contour(B.XC,B.YC,SMITH_mask,[0.5,0.5],'EdgeColor','k','LineWidth',2)
plot(SMITH_shape,'FaceColor','none','EdgeColor','k','LineWidth',2)
axis equal tight

% FUNCTIONS
function [region_mask, region_shape] = mask_region(polyvec_gr,polyvec_fl,Xq,Yq,mask)
region_shape_gr = union(polyvec_gr); % create union of the polyshape vector
region_mask_gr = inpolygon(Xq,Yq,region_shape_gr.Vertices(:,1),region_shape_gr.Vertices(:,2)) & mask; % mask the mesh center points for grounded nunataks

region_shape_fl = union(polyvec_fl); % create union of the polyshape vector
region_shape_fl = rmholes(region_shape_fl); % remove holes from polyshape
region_mask_fl = inpolygon(Xq,Yq,region_shape_fl.Vertices(:,1),region_shape_fl.Vertices(:,2)); % do not mask the mesh center points for floating ice

region_shape = union([region_shape_gr,region_shape_fl]); % create union of the polyshape vector
region_mask = (region_mask_gr | region_mask_fl); % mask of points in the shape
end

%%
t = all_requested_timesteps;
bed = md.geometry.bed;
for all_requested_timesteps
    thickness_matrix(:,i) = [md.results.TransientSolution(1).Thickness]; % add this time step to the matrix
end

thickness_dan = InterpFromMeshToMesh2d(md.mesh.elements, md.mesh.x, md.mesh.y, thickness_matrix, X, Y);