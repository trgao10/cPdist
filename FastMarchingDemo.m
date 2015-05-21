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



% path(pathdef);
% path(path, genpath('./utils'));
% 
% load('./samples/PNAS/Q13.mat');
% 
% %%% example of a trivial funtion
% f1 = ones(G.nV,1);
% G.ViewFunctionOnMesh(f1, []);
% 
% f2 = f1;
% f2(1:100) = -1;
% G.ViewFunctionOnMesh(f2, []);
% 
% close all
% % G.ViewFunctionOnMesh(G.Aux.VertArea',struct('mode','native'));
% % G.ViewFunctionOnMesh(G.Aux.Conf,struct('mode','native'));
% % localMaxInds = G.FindLocalMax(G.Aux.Conf,1,0);
% % hold on
% % scatter3(G.V(1,localMaxInds),G.V(2,localMaxInds),G.V(3,localMaxInds),...
% %     20,'k','filled');
% [Umin,Umax,Cmin,Cmax,Cmean,Cgauss,Normal] = G.ComputeCurvature();
% % G.ViewFunctionOnMesh(Cgauss,struct('mode','native'));
% % figure;G.ViewFunctionOnMesh(Cmean,struct('mode','native'));
% 
% G.draw();
% hold on
% quiver3(G.V(1,:),G.V(2,:),G.V(3,:),Normal(1,:),Normal(2,:),Normal(3,:));
% 
% H = Mesh('VF',G.V+0.1*Normal,G.F);
% H.draw();
