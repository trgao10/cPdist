%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% clean up meshes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% G = Mesh('off','../DATA/PNAS/meshes/j14_sas.off');
% NaNInds = 1780;
% 
% G.DeleteVertex(1780);
% 
% figure;G.draw();
% hold on;
% scatter3(G.V(1,NaNInds),G.V(2,NaNInds),G.V(3,NaNInds),30,'g','filled');
% % scatter3(G.V(1,NaNInds),G.V(2,NaNInds),G.V(3,NaNInds),30,'g','filled');
% G.Write('./j14_sas.off','off',[]);
% 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% compare distances and landmark MSEs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GroupSize = 116;
% 
% load('./results/Teeth/cPdist/cPdistMatrix.mat');
% figure;
% imagesc(cPdistMatrix./max(cPdistMatrix(:))*64);
% axis equal;
% axis([1,GroupSize,1,GroupSize]);
% 
% cPlmkMSEMatrix = load('./results/Teeth/cPdist/lmkMSEMatrix.mat');
% cPlmkMSEMatrix = cPlmkMSEMatrix.lmkMSEMatrix;
% figure;
% imagesc(cPlmkMSEMatrix./max(cPlmkMSEMatrix(:))*64);
% axis equal;
% axis([1,GroupSize,1,GroupSize]);
% 
% load('./results/Teeth/cPMST/FeatureFixOff/cPMSTDistMatrix.mat');
% figure;
% imagesc(ImprDistMatrix./max(cPdistMatrix(:))*64);
% axis equal;
% axis([1,GroupSize,1,GroupSize]);
% 
% cPMSTlmkMSEMatrix = load('./results/Teeth/cPMST/FeatureFixOff/cPMSTlmkMSEMatrix.mat');
% cPMSTlmkMSEMatrix = cPMSTlmkMSEMatrix.lmkMSEMatrix;
% figure;
% imagesc(cPMSTlmkMSEMatrix./max(cPMSTlmkMSEMatrix(:))*64);
% axis equal;
% axis([1,GroupSize,1,GroupSize]);
% 
% close all;
% 
% figure;
% scatter(cPdistMatrix(:),ImprDistMatrix(:),10,'g');
% axis equal;
% hold on;
% title('distances before/after MST improvement');
% plot([0,0.2],[0,0.2],'r');
% axis([0,0.2,0,0.2]);
% 
% figure;
% scatter(cPlmkMSEMatrix(:),cPMSTlmkMSEMatrix(:),10,'b');
% axis equal;
% hold on;
% title('landmark MSEs before/after MST improvement');
% plot([0,0.7],[0,0.7],'r');
% axis([0,0.7,0,0.7]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% test FeatureFix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% preparation
clear all;
close all;
path(pathdef);
addpath(path,genpath([pwd '/utils/']));

%%% setup paths
base_path = [pwd '/'];
data_path = '../DATA/PNAS/';
rslts_path = [base_path 'rslts/'];
cluster_path = [base_path 'cluster/'];
samples_path = [base_path 'samples/Teeth/'];
TextureCoords1_path = [pwd '/results/Teeth/cPMST/FeatureFixOff/TextureCoords1/'];
TextureCoords2_path = [pwd '/results/Teeth/cPMST/FeatureFixOff/TextureCoords2/'];
landmarks_path = [data_path 'landmarks_teeth.mat'];

% GroupSize = length(taxa_code);
chunk_size = 55;
rslt_mat = [rslts_path 'rslt_mat_' num2str(1)];

Names = {'D09','B03'};

G1 = [samples_path Names{1} '.mat'];
G2 = [samples_path Names{2} '.mat'];

options.TextureCoords1Path = TextureCoords1_path;
options.TextureCoords2Path = TextureCoords2_path;
options.ChunkSize = chunk_size;

GM = load(G1);
GM = GM.G;
GN = load(G2);
GN = GN.G;

load(rslt_mat);

tic;
disp(['Fixing features for ' GM.Aux.name ' vs ' GN.Aux.name '...']);
rslt = FeatureFix(GM,GN,2,1,options);
lk2 = GN.V(:,GetLandmarks(GN,LandmarksPath));
lk1 = GN.V(:,rslt.ImprMap(GetLandmarks(GM,LandmarksPath)));
rslt.lkMSE = mean(sqrt(sum((lk2-lk1).^2)));
save(rslt_mat,'Imprrslt');
disp(['Feature Fixing for ' GM.Aux.name ' vs ' GN.Aux.name ' done.']);
toc;






