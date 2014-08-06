%% preparation
clear vars;
close all;
path(pathdef);
addpath(path,genpath('./utils/'));

%% set parameters
obj_path = [pwd '/obj/'];
delete_command = 'rm -f ';
data_path = '~/Work/DATA/PNAS/';
meshes_path = [data_path 'meshes/'];
Names = {'w02','Q19'};

%% parse parameters
if ~exist(obj_path, 'dir')
    mkdir(obj_path);
end
command_text = [delete_command obj_path '*.obj'];
system(command_text);
disp(command_text);

taxa_code = load([data_path 'teeth_taxa_table.mat']);
taxa_code = taxa_code.taxa_code;

TAXAind1 = find(strcmpi(taxa_code, Names{1}));
TAXAind2 = find(strcmpi(taxa_code, Names{2}));

%% load mesh and uniformize
GM = Mesh('off',[meshes_path taxa_code{TAXAind1} '_sas.off']);
GN = Mesh('off',[meshes_path taxa_code{TAXAind2} '_sas.off']);

GM.Centralize('ScaleArea');
GN.Centralize('ScaleArea');

GM.ComputeMidEdgeUniformization;
GN.ComputeMidEdgeUniformization;

%% compute continuous Procrustes distance
[cPdist12,cPmap,TextureCoords1,TextureCoords2,ref12] = GM.ComputeContinuousProcrustes(GN);

%% print maps to texture coordinates
obj_surf_1 = [obj_path '1.obj'];
obj_surf_2 = [obj_path '2.obj'];
if ref12==1
    GM.F = GM.F([2,1,3],:);
end
options.Texture.Coordinates = TextureCoords1/2+0.5;
GM.Write(obj_surf_1,'obj',options);
options.Texture.Coordinates = TextureCoords2/2+0.5;
GN.Write(obj_surf_2,'obj',options);


