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
Names = {'a15','a16'};

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

options = [];

taxa_code = load([data_path 'teeth_taxa_table.mat']);
taxa_code = taxa_code.taxa_code;

TAXAind1 = find(strcmpi(taxa_code, Names{1}));
TAXAind2 = find(strcmpi(taxa_code, Names{2}));

%% load mesh and uniformize
if ~exist([sample_path Names{1} '.mat'],'file')
    GM = Mesh('off',[meshes_path taxa_code{TAXAind1} '_sas.off']);
    GM.Aux.name = Names{1};
    GM.Centralize('ScaleArea');
    GM.ComputeMidEdgeUniformization(options);
    G = Mesh(GM);
    save([sample_path Names{1} '.mat'], 'G');
else
    GM = load([sample_path Names{1} '.mat']);
    GM = GM.G;
end

if ~exist([sample_path Names{2} '.mat'],'file')
    GN = Mesh('off',[meshes_path taxa_code{TAXAind2} '_sas.off']);
    GN.Aux.name = Names{2};
    GN.Centralize('ScaleArea');
    GN.ComputeMidEdgeUniformization(options);
    G = Mesh(GN);
    save([sample_path Names{2} '.mat'], 'G');
else
    GN = load([sample_path Names{2} '.mat']);
    GN = GN.G;
end

% options.ConfMaxLocalWidth = 5;
% options.GaussMaxLocalWidth = 5;
% options.GaussMinLocalWidth = 5;
% options.ADMaxLocalWidth = 5;
% options.ExcludeBoundary = 1;
% options.Display = 'on';

%% compute continuous Procrustes distance
[cPdist12,cPmap12,TextureCoords1,TextureCoords2,ref12] = GM.ComputeContinuousProcrustes(GN);
[cPdist21,cPmap21,~,~,~] = GN.ComputeContinuousProcrustes(GM);

%% print maps to texture coordinates
obj_surf_1 = [obj_path '1.obj'];
obj_surf_2 = [obj_path '2.obj'];
options.Texture.Coordinates = TextureCoords1/2+0.5;
GM.Write(obj_surf_1,'obj',options);
options.Texture.Coordinates = TextureCoords2/2+0.5;
GN.Write(obj_surf_2,'obj',options);

%% visualize landmarks
[~,R,~] = MapToDist(GM.V,GN.V,cPmap12,GM.Aux.VertArea);
sGM = Mesh(GM);
sGM.V = R*GM.V;

options.type = 'full';
options.landmarks = 'on';
options.LandmarksPath = [data_path 'landmarks_teeth.mat'];
options.MeshesPath = [data_path 'meshes/'];
ViewTeethMapS(sGM, GN, {cPmap12,cPmap21}, options);

