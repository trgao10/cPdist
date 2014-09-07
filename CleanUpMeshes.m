%% Preparation
clear vars;
close all;
path(pathdef);
addpath(path,genpath([pwd '/utils/']));

%% set path
% meshesPath = '/media/trgao10/Work/MATLAB/DATA/PNAS/meshes/';
meshesPath = '/media/trgao10/Work/MATLAB/DATA/PNAS/newMeshes/';

%% read all meshes from meshesPath
meshFiles = dir(meshesPath);
meshFiles(1:2) = [];

%% set option for deleting isolated vertices
options.display = 'on';
options.exclude_boundary = 1;

for j=1:length(meshFiles)
    progressbar(j,length(meshFiles),20);
    G = Mesh('off',[meshesPath meshFiles(j).name]);
    [F,Inds] = unique(sort(G.F',2),'rows','first');
    if ~isempty(setdiff(1:G.nF,Inds))
        disp([meshFiles(j).name ' contains duplicate faces!']);
        V = G.V;
        F(setdiff(1:G.nF,Inds),:) = [];
        clear G;
        G = Mesh('VF',V,F');
    end
    dVInds = G.DeleteIsolatedVertex(options);
    if ~isempty(dVInds);
        disp([meshFiles(j).name ' contains non-boundary isolated vertex!']);
    end
%     G.Write([newMeshesPath meshFiles(j).name],'off',[]);
    clear G;
end

