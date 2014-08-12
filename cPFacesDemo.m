%% preparation
clear vars;
% close all;
path(pathdef);
addpath(path,genpath([pwd '/utils/']));

%% set parameters
Names = {'004','009'};
threshold = 0.1;

options.FeatureType = 'ConfMax';
options.NumDensityPnts = 100;
options.AngleIncrement = 0.05;
options.NumFeatureMatch = 4;
options.GaussMinMatch = 'off';

obj_path = [pwd '/obj/'];
sample_path = [pwd '/samples_face/'];
data_path = [pwd '/'];
meshes_path = data_path;
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


