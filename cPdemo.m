%% preparation
clear vars;
close all;
path(pathdef);
addpath(path,genpath('./utils/'));

%% set parameters
obj_path = [pwd '/obj/'];
sample_path = [pwd '/sample/'];
delete_command = 'rm -f ';
data_path = '~/Work/DATA/PNAS/';
meshes_path = [data_path 'meshes/'];
Names = {'Q18','Q19'};

%% parse parameters
if ~exist(obj_path, 'dir')
    mkdir(obj_path);
end
command_text = [delete_command obj_path '*.obj'];
system(command_text);
disp(command_text);
if ~exist(sample_path, 'dir')
    mkdir(sample_path);
end

Gs = cell(2,1);
options = [];

taxa_code = load([data_path 'teeth_taxa_table.mat']);
taxa_code = taxa_code.taxa_code;
TAXAind = cellfun(@(name) find(strcmpi(taxa_code,name)),Names);

%% load mesh and uniformize
for j=1:2
    if ~exist([sample_path Names{j} '.mat'],'file')
        Gs{j} = Mesh('off',[meshes_path taxa_code{TAXAind(j)} '_sas.off']);
        Gs{j}.Aux.name = Names{j};
        Gs{j}.Centralize('ScaleArea');
        Gs{j}.ComputeMidEdgeUniformization(options);
        G = Mesh(Gs{j});
        save([sample_path Names{j} '.mat'], 'G');
    else
        Gs{j} = load([sample_path Names{j} '.mat']);
        Gs{j} = Gs{j}.G;
    end
end

% options.ConfMaxLocalWidth = 5;
% options.GaussMaxLocalWidth = 5;
% options.GaussMinLocalWidth = 5;
% options.ADMaxLocalWidth = 5;
% options.ExcludeBoundary = 1;
% options.Display = 'on';

%% compute continuous Procrustes distance
rslt12 = Gs{1}.ComputeContinuousProcrustes(Gs{2});
rslt21 = Gs{2}.ComputeContinuousProcrustes(Gs{1});

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

%% visualize landmarks
[~,R,~] = MapToDist(Gs{1}.V,Gs{2}.V,rslt12.cPmap,Gs{1}.Aux.VertArea);
sGM = Mesh(Gs{1});
sGM.V = R*Gs{1}.V;

options.type = 'full';
options.landmarks = 'on';
options.LandmarksPath = [data_path 'landmarks_teeth.mat'];
options.MeshesPath = [data_path 'meshes/'];
ViewTeethMapS(sGM, Gs{2}, {rslt12.cPmap,rslt21.cPmap}, options);

