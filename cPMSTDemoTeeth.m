%% preparation
clear vars;
% close all;
path(pathdef);
addpath(path,genpath([pwd '/utils/']));

%% set parameters
Names = {'H16','j18'};

obj_path = [pwd '/obj/'];
sample_path = [pwd '/samples/Teeth/'];
cPmaps_path = [pwd '/results/Teeth/cPdist/cPmapsMatrix.mat'];
data_path = '~/Work/MATLAB/DATA/PNAS/';
delete_command = 'rm -f ';

%% parse parameters
touch(sample_path);
touch(obj_path);
command_text = [delete_command obj_path '*.obj'];
system(command_text);
disp(command_text);

Gs = cell(2,1);

taxa_code = load([data_path 'teeth_taxa_table.mat']);
taxa_code = taxa_code.taxa_code;
TAXAind = cellfun(@(name) find(strcmpi(taxa_code,name)),Names);

%% load all cPmaps (slow)
disp('loading all cPmaps...');
load(cPmaps_path);
disp('loaded');

%% load flattend meshes
for j=1:2
    Gs{j} = load([sample_path taxa_code{TAXAind(j)} '.mat']);
    Gs{j} = Gs{j}.G;
end

MSTpath = FindMSTPath(TAXAind{1},TAXAind{2},cPdistMatrix);
cPMSTmap12 = ComposeMapsAlongPath(MSTpath,maps_path);
cPMSTdist12 = MapToDist(Gs{1}.V,Gs{2}.V,cPMSTmap12,Gs{1}.Aux.VertArea);
cPMSTmap21 = ComposeMapsAlongPath(MSTpath(end:-1:1),maps_path);
cPMSTdist21 = MapToDist(Gs{2}.V,Gs{1}.V,cPMSTmap21,Gs{2}.Aux.VertArea);

%% print maps to texture coordinates
obj_surf_1 = [obj_path '1.obj'];
obj_surf_2 = [obj_path '2.obj'];
if cPMSTdist12<cPMSTdist21
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
[~,R,~] = MapToDist(Gs{1}.V,Gs{2}.V,cPMSTmap12,Gs{1}.Aux.VertArea);
sGM = Mesh(Gs{1});
sGM.V = R*Gs{1}.V;

options.type = 'full';
options.landmarks = 'on';
options.LandmarksPath = [data_path 'landmarks_teeth.mat'];
options.MeshesPath = [data_path 'meshes/'];
ViewTeethMapS(sGM, Gs{2}, {cPMSTmap12,cPMSTmap21}, options);

