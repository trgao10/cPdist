%% preparation
clear vars;
% close all;
path(pathdef);
addpath(path,genpath([pwd '/utils/']));

%% set parameters
Names = {'Angry/044','Happy/052'};
% threshold = 0.1;

options.FeatureType = 'ConfMax';
options.NumDensityPnts = 100;
options.AngleIncrement = 0.05;
options.NumFeatureMatch = 3;
options.GaussMinMatch = 'off';
options.SmoothCurvatureFields = 3;
options.DensityLocalWidth = 3;

options.ConfMaxLocalWidth = 3;
options.GaussMaxLocalWidth = 3;
options.GaussMinLocalWidth = 3;
options.ADMaxLocalWidth = 3;
options.ExcludeBoundary = 1;
options.Display = 'off';

obj_path = [pwd '/obj/'];
sample_path = [pwd '/samples/Faces/'];
data_path = '/media/trgao10/Work/MATLAB/Faces/';
meshes_path = [data_path 'F001/'];
delete_command = 'rm -f ';

%% parse parameters
touch(sample_path);
touch(obj_path);
command_text = [delete_command obj_path '*.obj'];
system(command_text);
disp(command_text);

Gs = cell(2,1);

%% load mesh and uniformize
for j=1:2
    if ~exist([sample_path Names{j} '.mat'],'file')
        Gs{j} = Mesh('off',[meshes_path Names{j} '.off']);
        Gs{j}.Aux.name = [meshes_path Names{j} '.off'];
        
        dVInds = Gs{j}.FindBoundaries;
        Gs{j}.DeleteVertex(dVInds);
        Gs{j}.DeleteIsolatedVertex;
        
%         dVInds = Gs{j}.FindBoundaries;
%         Gs{j}.DeleteVertex(dVInds);
%         Gs{j}.DeleteIsolatedVertex;
%         
%         dVInds = Gs{j}.FindBoundaries;
%         Gs{j}.DeleteVertex(dVInds);
%         Gs{j}.DeleteIsolatedVertex;
        
        Gs{j}.Centralize('ScaleArea');
        Gs{j}.ComputeMidEdgeUniformization(options);
        G = Mesh(Gs{j});
        save([sample_path Names{j} '.mat'], 'G');
    else
        Gs{j} = load([sample_path Names{j} '.mat']);
        Gs{j} = Gs{j}.G;
    end
end

%% compute continuous Procrustes distance
rslt12 = Gs{1}.ComputeContinuousProcrustes(Gs{2},options);
rslt21 = Gs{2}.ComputeContinuousProcrustes(Gs{1},options);

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


