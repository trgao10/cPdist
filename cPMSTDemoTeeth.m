%% Preparation
clear vars;
% close all;
path(pathdef);
addpath(path,genpath([pwd '/utils/']));

%% Set Parameters
Names = {'Q19','b08'};
% Names = {'H08','w01'};

options.ImprType = 'MST';
options.SmoothMap = 0;
options.Iter = 'on';

obj_path = [pwd '/obj/'];
sample_path = [pwd '/samples/Teeth/'];
cPmaps_path = [pwd '/results/Teeth/cPdist/cPmapsMatrix.mat'];
cPdist_path = [pwd '/results/Teeth/cPdist/cPdistMatrix.mat'];
data_path = '~/Work/MATLAB/DATA/PNAS/';
delete_command = 'rm -f ';

%%% options for ViewTeethMapS %%%
options.LandmarksPath = [data_path 'landmarks_teeth.mat'];
options.MeshesPath = [data_path 'meshes/'];

%% Parse Parameters
touch(sample_path);
touch(obj_path);
command_text = [delete_command obj_path '*.obj'];
system(command_text);
disp(command_text);

Gs = cell(2,1);

taxa_code = load([data_path 'teeth_taxa_table.mat']);
taxa_code = taxa_code.taxa_code;
TAXAind = cellfun(@(name) find(strcmpi(taxa_code,name)),Names);

%% Load All cPmaps and cPdists, View MST
load(cPdist_path);
% ConstructMST(cPdistMatrix,taxa_code);

disp('loading all cPmaps...');
load(cPmaps_path); % load cell arry "cPmapsMatrix"
disp('loaded');

%% Load Flattend Meshes
for j=1:2
    Gs{j} = load([sample_path taxa_code{TAXAind(j)} '.mat']);
    Gs{j} = Gs{j}.G;
end

%% Visualize Landmark Propagation for Original Maps
[~,R,~] = MapToDist(Gs{1}.V,Gs{2}.V,cPmapsMatrix{TAXAind(1),TAXAind(2)},Gs{1}.Aux.VertArea);
tGM = Mesh(Gs{1});
tGM.V = R*Gs{1}.V;

ViewTeethMapS(tGM, Gs{2}, {cPmapsMatrix{TAXAind(1),TAXAind(2)},cPmapsMatrix{TAXAind(2),TAXAind(1)}}, options);
set(gcf,'Name','cP');

%% Improve Maps
if strcmpi(options.ImprType,'MST')
    MSTpath = FindMSTPath(TAXAind(1),TAXAind(2),cPdistMatrix);
    disp('MST Path: ');
    disp(taxa_code(MSTpath));
end

rslt12 = Gs{1}.ImproveMap(Gs{2},cPdistMatrix,cPmapsMatrix,taxa_code,options);
rslt21 = Gs{2}.ImproveMap(Gs{1},cPdistMatrix,cPmapsMatrix,taxa_code,options);

%% Visualize Landmark Propagation for Improved Maps
[~,R,~] = MapToDist(Gs{1}.V,Gs{2}.V,rslt12.ImprMap,Gs{1}.Aux.VertArea);
sGM = Mesh(Gs{1});
sGM.V = R*Gs{1}.V;

ViewTeethMapS(sGM, Gs{2}, {rslt12.ImprMap,rslt21.ImprMap}, options);
set(gcf,'Name',options.ImprType);

%% Print Maps to Texture Coordinates
obj_surf_1 = [obj_path '1.obj'];
obj_surf_2 = [obj_path '2.obj'];

if rslt12.ImprDist<rslt21.ImprDist
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

% %% Pairly Improve Maps
% [rslt] = PairlyImproveMaps(Gs{1},Gs{2},rslt12,rslt21,options);
% 
% %% Visualize Landmark Propagation for Pairly Improved Maps
% [~,R,~] = MapToDist(Gs{1}.V,Gs{2}.V,rslt.ImprMap,Gs{1}.Aux.VertArea);
% pGM = Mesh(Gs{1});
% pGM.V = R*Gs{1}.V;
% 
% ViewTeethMapS(pGM, Gs{2}, {rslt.ImprMap,rslt.invImprMap}, options);
% set(gcf,'Name',[options.ImprType ' Pairly Improved']);
% 
% %% Print Maps to Texture Coordinates
% obj_surf_1 = [obj_path '1.obj'];
% obj_surf_2 = [obj_path '2.obj'];
% 
% options.Texture.Coordinates = rslt.TextureCoords1/2+0.5;
% Gs{1}.Write(obj_surf_1,'obj',options);
% options.Texture.Coordinates = rslt.TextureCoords2/2+0.5;
% Gs{2}.Write(obj_surf_2,'obj',options);

% if rslt12.ImprDist<rslt21.ImprDist
%     options.Texture.Coordinates = rslt12.TextureCoords1/2+0.5;
%     Gs{1}.Write(obj_surf_1,'obj',options);
%     options.Texture.Coordinates = rslt12.TextureCoords2/2+0.5;
%     Gs{2}.Write(obj_surf_2,'obj',options);
% else
%     options.Texture.Coordinates = rslt21.TextureCoords2/2+0.5;
%     Gs{1}.Write(obj_surf_1,'obj',options);
%     options.Texture.Coordinates = rslt21.TextureCoords1/2+0.5;
%     Gs{2}.Write(obj_surf_2,'obj',options);
% end

