%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% clean up meshes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% path(pathdef);
% addpath(path,genpath([pwd '/utils/']));
% 
% G = Mesh('off','../DATA/PNAS/meshes/s17_sas.off');
% NaNInds = 2466;
% 
% G.DeleteVertex(NaNInds);
% 
% figure;G.draw();
% hold on;
% scatter3(G.V(1,NaNInds),G.V(2,NaNInds),G.V(3,NaNInds),30,'g','filled');
% G.Write('./s17_sas.off','off',[]);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% compare distances and landmark MSEs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
GroupSize = 116;

load('./results/Teeth/cPDist/cPDistMatrix.mat');
% figure;
% imagesc(cPDistMatrix./max(cPDistMatrix(:))*64);
% axis equal;
% axis([1,GroupSize,1,GroupSize]);

load('./results/Teeth/cPDist/cPlmkMSEMatrix.mat');
% figure;
% imagesc(cPlmkMSEMatrix./max(cPlmkMSEMatrix(:))*64);
% axis equal;
% axis([1,GroupSize,1,GroupSize]);

cPMSTFeatureFixOff = load('./results/Teeth/cPMST/FeatureFixOff/cPMSTDistMatrix.mat');
cPMSTFeatureFixOff = cPMSTFeatureFixOff.ImprDistMatrix;
cPMSTFFofflmkMSEMatrix = load('./results/Teeth/cPMST/FeatureFixOff/cPMSTlmkMSEMatrix.mat');
cPMSTFFofflmkMSEMatrix = cPMSTFFofflmkMSEMatrix.lmkMSEMatrix;
cPMSTFeatureFixOn = load('./results/Teeth/cPMST/FeatureFixOn/cPMSTDistMatrix.mat');
cPMSTFeatureFixOn = cPMSTFeatureFixOn.ImprDistMatrix;
cPMSTFFonlmkMSEMatrix = load('./results/Teeth/cPMST/FeatureFixOn/cPMSTlmkMSEMatrix.mat');
cPMSTFFonlmkMSEMatrix = cPMSTFFonlmkMSEMatrix.lmkMSEMatrix;

cPViterbiFeatureFixOff = load('./results/Teeth/cPViterbi/FeatureFixOff/cPViterbiDistMatrix.mat');
cPViterbiFeatureFixOff = cPViterbiFeatureFixOff.ImprDistMatrix;
cPViterbiFFofflmkMSEMatrix = load('./results/Teeth/cPViterbi/FeatureFixOff/cPViterbilmkMSEMatrix.mat');
cPViterbiFFofflmkMSEMatrix = cPViterbiFFofflmkMSEMatrix.lmkMSEMatrix;
cPViterbiFeatureFixOn = load('./results/Teeth/cPViterbi/FeatureFixOn/cPViterbiDistMatrix.mat');
cPViterbiFeatureFixOn = cPViterbiFeatureFixOn.ImprDistMatrix;
cPViterbiFFonlmkMSEMatrix = load('./results/Teeth/cPViterbi/FeatureFixOn/cPViterbilmkMSEMatrix.mat');
cPViterbiFFonlmkMSEMatrix = cPViterbiFFonlmkMSEMatrix.lmkMSEMatrix;

cPLASTFeatureFixOff = load('./results/Teeth/cPLAST/FeatureFixOff/cPLASTDistMatrix.mat');
cPLASTFeatureFixOff = cPLASTFeatureFixOff.ImprDistMatrix;
cPLASTFFofflmkMSEMatrix = load('./results/Teeth/cPLAST/FeatureFixOff/cPLASTlmkMSEMatrix.mat');
cPLASTFFofflmkMSEMatrix = cPLASTFFofflmkMSEMatrix.lmkMSEMatrix;
cPLASTFeatureFixOn = load('./results/Teeth/cPLAST/FeatureFixOn/cPLASTDistMatrix.mat');
cPLASTFeatureFixOn = cPLASTFeatureFixOn.ImprDistMatrix;
cPLASTFFonlmkMSEMatrix = load('./results/Teeth/cPLAST/FeatureFixOn/cPLASTlmkMSEMatrix.mat');
cPLASTFFonlmkMSEMatrix = cPLASTFFonlmkMSEMatrix.lmkMSEMatrix;

figure;
scatter(cPLASTFeatureFixOff(:),cPLASTFeatureFixOn(:),10,'g');
axis equal;
hold on;
title('cPViterbi distances before/after FeatureFix');
plot([0,0.2],[0,0.2],'r');
axis([0,0.2,0,0.2]);

figure;
scatter(cPLASTFFofflmkMSEMatrix(:),cPLASTFFonlmkMSEMatrix(:),10,'b');
axis equal;
hold on;
title('cPViterbi landmark MSEs before/after FeatureFix');
plot([0,0.7],[0,0.7],'r');
axis([0,0.7,0,0.7]);

% close all;
% figure;
% scatter(cPDistMatrix(:),ImprDistMatrix(:),10,'g');
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
% %%% preparation
% clear all;
% close all;
% path(pathdef);
% addpath(path,genpath([pwd '/utils/']));
% 
% %%% setup paths
% base_path = [pwd '/'];
% data_path = '../DATA/PNAS/';
% rslts_path = [base_path 'rslts/'];
% cluster_path = [base_path 'cluster/'];
% samples_path = [base_path 'samples/Teeth/'];
% TextureCoords1_path = [pwd '/results/Teeth/cPViterbi/FeatureFixOff/TextureCoords1/'];
% TextureCoords2_path = [pwd '/results/Teeth/cPViterbi/FeatureFixOff/TextureCoords2/'];
% LandmarksPath = [data_path 'landmarks_teeth.mat'];
% TaxaCode_path = [data_path 'teeth_taxa_table.mat'];
% 
% %%% load taxa codes
% taxa_code = load(TaxaCode_path);
% TaxaCode = taxa_code.taxa_code;
% GroupSize = length(taxa_code);
% ChunkSize = 55;
% 
% Names = {'s17','J12'};
% 
% G1 = [samples_path Names{1} '.mat'];
% G2 = [samples_path Names{2} '.mat'];
% 
% options.TextureCoords1Path = TextureCoords1_path;
% options.TextureCoords2Path = TextureCoords2_path;
% options.ChunkSize = ChunkSize;
% 
% GM = load(G1);
% GM = GM.G;
% GN = load(G2);
% GN = GN.G;
% 
% TAXAind = cellfun(@(name) find(strcmpi(TaxaCode,name)),{GM.Aux.name,GN.Aux.name});
% TAXAind1 = TAXAind(1);
% TAXAind2 = TAXAind(2);
% 
% ChunkIdx = @(TAXAind1,TAXAind2) ceil(((TAXAind1-1)*GroupSize+TAXAind2)/ChunkSize);
% rslt_mat = [rslts_path 'rslt_mat_' num2str(ChunkIdx(TAXAind1,TAXAind2))];
% load(rslt_mat);
% 
% tic;
% disp(['Fixing features for ' GM.Aux.name ' vs ' GN.Aux.name '...']);
% rslt = FeatureFix(GM,GN,TAXAind1,TAXAind2,options);
% lk2 = GN.V(:,GetLandmarks(GN,LandmarksPath));
% lk1 = GN.V(:,rslt.ImprMap(GetLandmarks(GM,LandmarksPath)));
% rslt.lkMSE = mean(sqrt(sum((lk2-lk1).^2)));
% % save(rslt_mat,'Imprrslt');
% disp(['Feature Fixing for ' GM.Aux.name ' vs ' GN.Aux.name ' done.']);
% toc;


