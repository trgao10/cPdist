%% preparation
clear vars;
% close all;
path(pathdef);
addpath(path,genpath([pwd '/utils/']));

%% set parameters
Names = {'MN-Rio-8513_M1220','AMNH-M-76882_M786'};

options.FeatureType = 'ConfMax';
options.NumDensityPnts = 1000;
options.AngleIncrement = 0.01;
% options.NumDensityPnts = 100;
% options.AngleIncrement = 0.05;
options.NumFeatureMatch = 4;
options.GaussMinMatch = 'off';
% options.FeatureMatchType = 'Laplacian';
% options.ConfMaxLocalWidth = 5;
% options.GaussMaxLocalWidth = 5;
% options.GaussMinLocalWidth = 5;
% options.ADMaxLocalWidth = 5;
% options.ExcludeBoundary = 1;
% options.Display = 'on';

obj_path = [pwd '/obj/'];
sample_path = [pwd '/samples/HDM/'];
data_path = '~/Work/MATLAB/DATA/HDM/';
meshes_path = [data_path 'meshes/'];
delete_command = 'rm -f ';

%% parse parameters
touch(sample_path);
touch(obj_path);
delete([obj_path '*.obj']);
% command_text = [delete_command obj_path '*.obj'];
% system(command_text);
% disp(command_text);

Gs = cell(2,1);

%% load mesh and uniformize
for j=1:2
    if ~exist([sample_path Names{j} '.mat'],'file')
        Gs{j} = Mesh('off',[meshes_path Names{j} '.off']);
        Gs{j}.Aux.name = Names{j};
        [Gs{j}.Aux.Area,Gs{j}.Aux.Center] = Gs{j}.Centralize('ScaleArea');
        Gs{j}.ComputeMidEdgeUniformization(options);
        Gs{j}.Nf = Gs{j}.ComputeFaceNormals;
        Gs{j}.Nv = Gs{j}.F2V'*Gs{j}.Nf';
        Gs{j}.Nv = Gs{j}.Nv'*diag(1./sqrt(sum((Gs{j}.Nv').^2,1)));
        Gs{j}.Aux.LB = Gs{j}.ComputeCotanLaplacian;
        G = Mesh(Gs{j});
        save([sample_path Names{j} '.mat'], 'G');
    else
        Gs{j} = load([sample_path Names{j} '.mat']);
        Gs{j} = Gs{j}.G;
    end
end

%% compute continuous Procrustes distance
tic;rslt12 = Gs{1}.ComputeContinuousProcrustes(Gs{2},options);toc;
tic;rslt21 = Gs{2}.ComputeContinuousProcrustes(Gs{1},options);toc;

disp(['rslt12.cPdist = ' num2str(rslt12.cPdist)]);
disp(['rslt21.cPdist = ' num2str(rslt21.cPdist)]);

%% print maps to texture coordinates
obj_surf_1 = [obj_path '1.obj'];
obj_surf_2 = [obj_path '2.obj'];
if rslt12.cPdist<rslt21.cPdist
    options.Texture.Coordinates = rslt12.TextureCoords1/2+0.5;
    Gs{1}.Write(obj_surf_1,'obj',options);
    options.Texture.Coordinates = rslt12.TextureCoords2/2+0.5;
    Gs{2}.Write(obj_surf_2,'obj',options);
else
    options.Texture.Coordinates = rslt21.TextureCoords2/2+0.5;
    Gs{1}.Write(obj_surf_1,'obj',options);
    options.Texture.Coordinates = rslt21.TextureCoords1/2+0.5;
    Gs{2}.Write(obj_surf_2,'obj',options);
end

