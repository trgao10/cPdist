clearvars;
close all;
path(pathdef);
path(path, genpath('./utils'));

samples_path = './samples/PNAS/';
MeshName = 'k11';
load([samples_path MeshName '.mat']);
% G = Mesh('off','./meshes/david-head.off');

%%% pick vertices of interest
%%% by invoking the picking function of the display
% G.draw();
% keyboard

%%% fast marching from initial points of interest
% D = G.PerformFastMarching(1707);
% G.ViewFunctionOnMesh(D,struct('mode','native'));

% Subsamples = G.GeodesicFarthestPointSampling(400);

%%% geodesic extraction demo
% geodesic = G.ComputeGeodesic(D, 1504);
% hold on
% plot3(geodesic(1,:),geodesic(2,:),geodesic(3,:),'k-','LineWidth',2);

% %%% anisotropic fast marching
% [K,M] = G.ComputeCurvature;
% % G.ViewFunctionOnMesh(M,struct('mode','native'));
% D = G.PerformFastMarching(1707,struct('W',M-min(M)));
% G.ViewFunctionOnMesh(D,struct('mode','native'));

% close all;
% KM = G.ComputeCPMS();
% DM = KM.ComputeUniformization(struct('method','LSCM','boundary_conditions','disc'));

NumPts = 100;
rng(100)
samples = zeros(1,NumPts);

samples(1) = randi(G.nV);

cback = 0;
for j=2:NumPts
    fprintf(repmat('\b',1,cback));
    cback = fprintf(['%2d/' num2str(NumPts) '\n'], j);
%     progressbar(j,NumPts);
    %%%% get samples(j) via FastMarching
    [~,samples(j)] = max(G.PerformFastMarching(samples(1:(j-1))));
end

[D,S,Q] = G.PerformFastMarching(samples);

G.draw();
hold on
scatter3(G.V(1,samples),G.V(2,samples),G.V(3,samples),20,'g','filled');

%%%%%% compute mesh parametrization
KM = G.ComputeCPMS(); %% curvature prescription mesh scaling
DM = KM.ComputeUniformization(struct('method','LSCM','boundary_conditions','disc'));

