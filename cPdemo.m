%% preparation
clear vars;
% close all;
path(pathdef);
addpath(path,genpath([pwd '/utils/']));

%% set parameters
Names = {'i14','b19'};
options.FeatureType = 'ConfMax';
options.NumDensityPnts = 100;
options.AngleIncrement = 0.05;
options.NumFeatureMatch = 4;
options.GaussMinMatch = 'off';
% options.ConfMaxLocalWidth = 5;
% options.GaussMaxLocalWidth = 5;
% options.GaussMinLocalWidth = 5;
% options.ADMaxLocalWidth = 5;
% options.ExcludeBoundary = 1;
% options.Display = 'on';

obj_path = [pwd '/obj/'];
sample_path = [pwd '/sample/'];
data_path = '~/Work/MATLAB/DATA/PNAS/';
meshes_path = [data_path 'meshes/'];
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

%% load mesh and uniformize
for j=1:2
    if ~exist([sample_path taxa_code{TAXAind(j)} '.mat'],'file')
        Gs{j} = Mesh('off',[meshes_path taxa_code{TAXAind(j)} '_sas.off']);
        Gs{j}.Aux.name = taxa_code{TAXAind(j)};
        Gs{j}.Centralize('ScaleArea');
        Gs{j}.ComputeMidEdgeUniformization(options);
        G = Mesh(Gs{j});
        save([sample_path taxa_code{TAXAind(j)} '.mat'], 'G');
    else
        Gs{j} = load([sample_path taxa_code{TAXAind(j)} '.mat']);
        Gs{j} = Gs{j}.G;
    end
end

%% compute continuous Procrustes distance
tic;
rslt12 = Gs{1}.ComputeContinuousProcrustes(Gs{2},options);
toc;
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

%% visualize landmarks
[~,R,~] = MapToDist(Gs{1}.V,Gs{2}.V,rslt12.cPmap,Gs{1}.Aux.VertArea);
sGM = Mesh(Gs{1});
sGM.V = R*Gs{1}.V;

options.type = 'full';
options.landmarks = 'on';
options.LandmarksPath = [data_path 'landmarks_teeth.mat'];
options.MeshesPath = [data_path 'meshes/'];
ViewTeethMapS(sGM, Gs{2}, {rslt12.cPmap,rslt21.cPmap}, options);

% %% iterative update
% [Dist12,proj_map12,VM12,VN12] = IterProjMoebius(Gs{1},Gs{2},rslt12.cPmap,rslt12.ref,options);
% [Dist21,proj_map21,VN21,VM21] = IterProjMoebius(Gs{2},Gs{1},rslt21.cPmap,rslt21.ref,options);
% 
% while Dist12~=Dist21
%     if Dist12<Dist21
%         inv_map12 = knnsearch(VM12',VN12');
%         [Dist21,proj_map21,VN21,VM21] = IterProjMoebius(Gs{2},Gs{1},inv_map12,rslt21.ref,options);
%         if Dist12<Dist21
%             break;
%         end
%     elseif Dist12>Dist21
%         inv_map21 = knnsearch(VN21',VM21');
%         [Dist12,proj_map12,VM12,VN12] = IterProjMoebius(Gs{1},G{2},inv_map21,rslt12.ref,options);
%         if Dist12>Dist21
%             break;
%         end
%     else
%         break;
%     end
% end
% 
% options.GaussMinInds = 'on';
% while Dist12~=Dist21
%     if Dist12<Dist21
%         inv_map12 = knnsearch(VM12',VN12');
%         [Dist21,proj_map21,VN21,VM21] = IterProjMoebius(Gs{2},Gs{1},inv_map12,rslt21.ref,options);
%         if Dist12<Dist21
%             break;
%         end
%     elseif Dist12>Dist21
%         inv_map21 = knnsearch(VN21',VM21');
%         [Dist12,proj_map12,VM12,VN12] = IterProjMoebius(Gs{1},Gs{2},inv_map21,rslt12.ref,options);
%         if Dist12>Dist21
%             break;
%         end
%     else
%         break;
%     end
% end
% 
% %++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% % summary and visualization
% %++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% disp(['Final Dist(' Names{1} ',' Names{2} ') = ' num2str(Dist12)]);
% disp(['Final Dist(' Names{2} ',' Names{1} ') = ' num2str(Dist21)]);
% maps = {proj_map12,proj_map21};
% 
% [~,R,~] = MapToDist(Gs{1}.V,Gs{2}.V,proj_map12,Gs{1}.Aux.VertArea);
% sGM = Mesh(Gs{1});
% sGM.V = R*Gs{1}.V;
% 
% options.type = 'full';
% options.landmarks = 'on';
% options.LandmarksPath = [data_path 'landmarks_teeth.mat'];
% options.MeshesPath = [data_path 'meshes/'];
% ViewTeethMapS(sGM, Gs{2}, maps, options);
% 
