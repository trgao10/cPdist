clearvars;
close all;
path(pathdef);
path(path, genpath('./utils'));

Names = {'U02','U05'};

G1 = load(['./samples/PNAS/' Names{1} '.mat']); G1 = G1.G;
G2 = load(['./samples/PNAS/' Names{2} '.mat']); G2 = G2.G;

figure;G1.draw();
figure;G2.draw();

%%%%%% compute mesh parametrization
K1 = G1.ComputeCPMS(); %% curvature prescription mesh scaling
D1 = K1.ComputeUniformization(struct('method','LSCM','boundary_conditions','disc'));

K2 = G2.ComputeCPMS(); %% curvature prescription mesh scaling
D2 = K2.ComputeUniformization(struct('method','LSCM','boundary_conditions','disc'));

G1.Aux.UniformizationV = D1.V;
G2.Aux.UniformizationV = D2.V;

options = struct('FeatureType','ConfMax','NumDensityPnts',100,...
    'AngleIncrement',0.05,'NumFeatureMatch',4);
rslt12 = G1.ComputeContinuousProcrustes(G2,options);

keyboard

%%%%%% write texture coordinates into .obj files
options.Texture.Coordinates = rslt12.TextureCoords1(1:2,:)/2+0.5;
G1.Write([Names{1} '.obj'],'obj',options);
options.Texture.Coordinates = rslt12.TextureCoords2(1:2,:)/2+0.5;
G2.Write([Names{2} '.obj'],'obj',options);
