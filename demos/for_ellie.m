%% preparation
clear vars;
% close all;
path(pathdef);
addpath(path,genpath('~/Work/MATLAB/cPdist/utils/'));

options = struct('FeatureType', 'ConfMax', 'NumDensityPnts', 100,...
    'AngleIncrement', 0.05, 'NumFeatureMatch', 4, 'GaussMinMatch', 'on');

load('../from_ellie/Sai_gc.mat');

for j=1:10
    ConfMax = Sai{j}.V(:,Sai{j}.Aux.ConfMaxInds);
    GaussMax = Sai{j}.V(:,Sai{j}.Aux.GaussMaxInds);
    GaussMin = Sai{j}.V(:,Sai{j}.Aux.GaussMinInds);
    ADMax = Sai{j}.V(:,Sai{j}.Aux.ADMaxInds);
    
    dVInds = Sai{j}.DeleteIsolatedVertex();
    if ~isempty(dVInds) %% recompute uniformization
        [Sai{j}.Aux.Area,Sai{j}.Aux.Center] = Sai{j}.Centralize('ScaleArea');
        ConfMax = (ConfMax-repmat(Sai{j}.Aux.Center,1,size(ConfMax,2)))/sqrt(Sai{j}.Aux.Area);
        GaussMax = (GaussMax-repmat(Sai{j}.Aux.Center,1,size(GaussMax,2)))/sqrt(Sai{j}.Aux.Area);
        GaussMin = (GaussMin-repmat(Sai{j}.Aux.Center,1,size(GaussMin,2)))/sqrt(Sai{j}.Aux.Area);
        ADMax = (ADMax-repmat(Sai{j}.Aux.Center,1,size(ADMax,2)))/sqrt(Sai{j}.Aux.Area);
        
        Sai{j}.ComputeMidEdgeUniformization(options);
        Sai{j}.Nf = Sai{j}.ComputeFaceNormals;
        Sai{j}.Nv = Sai{j}.F2V'*Sai{j}.Nf';
        Sai{j}.Nv = Sai{j}.Nv'*diag(1./sqrt(sum((Sai{j}.Nv').^2,1)));
        Sai{j}.Aux.LB = Sai{j}.ComputeCotanLaplacian;
        
        TREE = kdtree_build(Sai{j}.V');
        Sai{j}.Aux.ConfMaxInds = kdtree_nearest_neighbor(TREE, ConfMax');
        Sai{j}.Aux.GaussMaxInds = kdtree_nearest_neighbor(TREE, GaussMax');
        Sai{j}.Aux.GaussMinInds = kdtree_nearest_neighbor(TREE, GaussMin');
        Sai{j}.Aux.ADMaxInds = kdtree_nearest_neighbor(TREE, ADMax');
    end
end

CPDMat = zeros(10);
CPRsltMat = cell(10);
for j=1:10
    for k=1:10
        CPRsltMat{j,k} = Sai{j}.ComputeContinuousProcrustes(Sai{k},options);
        CPDMat(j,k) = CPRsltMat{j,k}.cPdist;
    end
end

% options.Texture.Coordinates = CPRsltMat{1,10}.TextureCoords1/2+0.5;
% Sai{1}.Write('../1.obj','obj',options);
% options.Texture.Coordinates = CPRsltMat{1,10}.TextureCoords2/2+0.5;
% Sai{10}.Write('../2.obj','obj',options);

imagesc(min(CPDMat,CPDMat'));
axis equal
axis([0.5,10.5,0.5,10.5])
colorbar


