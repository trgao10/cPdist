%% preparation
clear vars;
% close all;
path(pathdef);
addpath(path,genpath('~/Work/MATLAB/cPdist/utils/'));

options = struct('FeatureType', 'GaussMax', 'NumDensityPnts', 100,...
    'AngleIncrement', 0.05, 'NumFeatureMatch', 4, 'GaussMinMatch', 'on');

GM = load('../samples/PNAS/ah16.mat');
GM = GM.G;

if ~exist('../k01_j18_05.mat', 'file')
    GN = Mesh('off', '../k01_j18_05.off');
    
    [GN.Aux.Area,GN.Aux.Center] = GN.Centralize('ScaleArea');
    GN.ComputeMidEdgeUniformization(options);
    GN.Nf = GN.ComputeFaceNormals;
    GN.Nv = GN.F2V'*GN.Nf';
    GN.Nv = GN.Nv'*diag(1./sqrt(sum((GN.Nv').^2,1)));
    GN.Aux.LB = GN.ComputeCotanLaplacian;
    G = Mesh(GN);
    save('../k01_j18_05.mat', 'G');
else
    GN = load('../k01_j18_05.mat');
    GN = GN.G;
end

rslt = GM.ComputeContinuousProcrustes(GN, options);

options.Texture.Coordinates = rslt.TextureCoords1/2+0.5;
GM.Write('../1.obj','obj',options);
options.Texture.Coordinates = rslt.TextureCoords2/2+0.5;
GN.Write('../2.obj','obj',options);



