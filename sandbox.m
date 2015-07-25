%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% compare two distance matrices (Lem/Gs)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('/home/trgao10/Work/MATLAB/cPdist/from_ellie/LemRslt_midedge.mat');
load('/home/trgao10/Work/MATLAB/cPdist/from_ellie/GsRslt_midedge.mat');

D1 = min(GsDistMat,GsDistMat'); %% symmetrize the distance matrix
D1 = D1-diag(diag(D1)); %% set diagonal elements to 0, avoid artifacts
d1 = squareform(D1, 'tovector'); %% extract the off-diagonal elements

D2 = min(LemDistMat,LemDistMat');
D2 = D2-diag(diag(D2));
d2 = squareform(D2, 'tovector');

hist([d1',d2']);

[h,p,ci,stats] = ttest2(d1,d2); %% two-sample t-test

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% compute cP distance from Ellie's downsampled meshes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% path(pathdef);
% addpath(path,genpath([pwd '/utils/']));
% 
% load('from_ellie/Lem.mat');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% view PCA scores as coordinates
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% result_path = '/';
% Roots = {'RootB03', 'RootMinSqDist'};
% Methods = {'cPComposedLAST', 'cPComposedLASTbalance', 'cPComposedLASTmean',...
%     'cPComposedLASTmedian', 'cPLAST', 'cPLASTbalance', 'cPMST', 'cPViterbi',...
%     'cPViterbiAngle0.5', 'cPViterbiAngle0.25'};
% FeatureFix = {'FeatureFixOff', 'FeatureFixOn'};
% NumPts = {'64','256','1024'};
% colors = ['r', 'g', 'b', 'm', 'k', 'y'];
% 
% RootsIdx = 1;
% MethodsIdx = 5;
% FeatureFixIdx = 1;
% 
% figure;
% % hold on
% for j=1:length(NumPts)
%     filePath = [result_path Roots{RootsIdx} filesep Methods{MethodsIdx} filesep FeatureFix{FeatureFixIdx} filesep];
%     if exist([filePath filesep 'RESULTS_' Methods{MethodsIdx} '_' FeatureFix{FeatureFixIdx} '_morphologika_' NumPts{j} '.csv'], 'file')
%         resultFileName = ['RESULTS_' Methods{MethodsIdx} '_' FeatureFix{FeatureFixIdx} '_morphologika_' NumPts{j} '.csv'];
%     else
%         resultFileName = [Methods{MethodsIdx} '_' FeatureFix{FeatureFixIdx} '_morphologika_' NumPts{j} '_RESULTS' '.csv'];
%     end
%     disp(resultFileName);
%     M = csvread([filePath resultFileName],687,2,[687,2,802,4]);
%     scatter3(M(:,1),M(:,2),M(:,3),20,colors(j),'filled');
%     hold on;
% end
% 
% M = csvread('./PNAS_Observer_Landmarks_18_RESULTS.csv',481,2,[481,2,596,4]);
% scatter3(M(:,1),M(:,2),M(:,3),20,'k','filled');
% 
% legend('64 pts', '256 pts', '1024 pts', 'Obs 18 pts');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% processing morphologika analysis results and write to morphologika file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% result_path = '/media/trgao10/Downloads/BoyerLab/';
% [~,TaxaCodes] = xlsread([result_path 'mergedResults.xlsx'], 'A03:A123');
% coordinates = xlsread([result_path 'mergedResults.xlsx'], 'B03:MK123');
% 
% %% Write Information to Output File
% fid = fopen('mergedResults_PC12_morphologika.txt','wt');
% if(fid==-1)
%     error('Can''t open the file.');
% end
% 
% %%% header
% fprintf(fid, '[Individuals]\n');
% fprintf(fid, '%d\n', length(TaxaCodes));
% fprintf(fid, '[landmarks]\n');
% fprintf(fid, '%d\n', size(coordinates,2)/3);
% fprintf(fid, '[dimensions]\n');
% fprintf(fid, '3\n');
% fprintf(fid, '[names]\n');
% for j=1:length(TaxaCodes)
%     fprintf(fid, '%s\n', TaxaCodes{j});
% end
% 
% fprintf(fid, '\n[rawpoints]\n');
% for j=1:length(TaxaCodes)
%     fprintf(fid,'\n%s\n\n', ['''' TaxaCodes{j}]);
%     reshapedCoords = reshape(coordinates(j,:),3,size(coordinates,2)/3);
%     reshapedCoords(3,:) = 0;
%     fprintf(fid, '%d %d %d\n', reshapedCoords);
% end
% 
% fclose(fid);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% processing morphologika analysis results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% result_path = '/media/trgao10/Downloads/BoyerLab/';
% Roots = {'RootB03', 'RootMinSqDist'};
% Methods = {'cPComposedLAST', 'cPComposedLASTbalance', 'cPComposedLASTmean',...
%     'cPComposedLASTmedian', 'cPLAST', 'cPLASTbalance', 'cPMST', 'cPViterbi',...
%     'cPViterbiAngle0.5', 'cPViterbiAngle0.25'};
% FeatureFix = {'FeatureFixOff', 'FeatureFixOn'};
% NumPts = {'64','256','1024'};
% 
% resultFileName = ['RESULTS_' Methods{1} '_' FeatureFix{1} '_morphologika_' NumPts{1} '.csv'];
% [~,TaxaCodes] = xlsread([result_path Roots{1} filesep Methods{1} filesep FeatureFix{1} filesep resultFileName],'B688:B803');
% TaxaCodes = (cellfun(@(x) strtrim(x), TaxaCodes, 'UniformOutput', false));
% Header = [TaxaCodes'; cell(2,length(TaxaCodes))];
% Header = reshape(Header, 1, numel(Header));
% Header = [Header; repmat({'PC1','PC2','PC3'},1,length(TaxaCodes))];
% Header = [cell(2,1), Header];
% 
% for j1=1:length(Roots)
%     for j2=1:length(Methods)
%         for j3=1:length(FeatureFix)
%             for j4=1:length(NumPts)
%                 filePath = [result_path Roots{j1} filesep Methods{j2} filesep FeatureFix{j3} filesep];
%                 if exist([filePath filesep 'RESULTS_' Methods{j2} '_' FeatureFix{j3} '_morphologika_' NumPts{j4} '.csv'], 'file')
%                     resultFileName = ['RESULTS_' Methods{j2} '_' FeatureFix{j3} '_morphologika_' NumPts{j4} '.csv'];
%                 else
%                     resultFileName = [Methods{j2} '_' FeatureFix{j3} '_morphologika_' NumPts{j4} '_RESULTS' '.csv'];
%                 end
%                 disp(resultFileName);                
%                 M = csvread([filePath resultFileName],687,2,[687,2,802,4]);
%                 Header = [Header;[[Methods{j2} '_' FeatureFix{j3} '_morphologika_' NumPts{j4} '_' Roots{j1}], num2cell(reshape(M',1,numel(M)))]];
%             end
%         end
%     end
% end
% 
% xlswrite('mergedResults.xlsx', Header);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% compare & contrast cPdistances/landmark MSEs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% result_path = '/media/trgao10/Work/MATLAB/ArchivedResults/Teeth/';
% % cPDistPath = [result_path 'cPDist/cPlmkMSEMatrix.mat'];
% % cPMSTDistOffPath = [result_path 'cPMST/FeatureFixOff/cPMSTlmkMSEMatrix.mat'];
% % cPMSTDistOnPath = [result_path 'cPMST/FeatureFixOn/cPMSTlmkMSEMatrix.mat'];
% cPDistPath = [result_path 'cPDist/cPDistMatrix.mat'];
% cPMSTDistOffPath = [result_path 'cPMST/FeatureFixOff/cPMSTDistMatrix.mat'];
% cPMSTDistOnPath = [result_path 'cPMST/FeatureFixOn/cPMSTDistMatrix.mat'];
% GroupSize = 116;
% 
% A = load(cPDistPath);
% A = A.cPDistMatrix;
% B = load(cPMSTDistOffPath);
% B = B.ImprDistMatrix;
% C = load(cPMSTDistOnPath);
% C = C.ImprDistMatrix;
% 
% figure;
% subplot(1,3,1);
% imagesc(A./max(A(:))*64);
% axis equal;
% axis([1,GroupSize,1,GroupSize]);
% title('cPDist');
% subplot(1,3,2);
% imagesc(B./max(B(:))*64);
% axis equal;
% axis([1,GroupSize,1,GroupSize]);
% title('cPMST Feature Fix Off');
% subplot(1,3,3);
% imagesc(C./max(C(:))*64);
% axis equal;
% axis([1,GroupSize,1,GroupSize]);
% title('cPMST Feature Fix On');
% 
% DistanceMatrices = {A,B,C};
% TitleTexts = {'cPDist','cPMST Feature Fix Off', 'cPMST Feature Fix On'};
% figure;
% count = 0;
% for i=1:3
%     for j=i+1:3
%         count = count+1;
%         subplot(1,3,count);
%         plot(DistanceMatrices{i}(:),DistanceMatrices{j}(:),'b.');
%         axis equal;
%         hold on;
%         axis([0,max(DistanceMatrices{i}(:)),0,max(DistanceMatrices{j}(:))]);
%         minLabel = min(max(DistanceMatrices{i}(:)),max(DistanceMatrices{j}(:)));
%         plot([0,minLabel],[0,minLabel],'r-');
%         title([TitleTexts{i} ' vs ' TitleTexts{j}]);
%     end
% end

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
% GroupSize = 116;
% ResultPath = '/media/trgao10/Work/MATLAB/ArchivedResults/Teeth/';
% 
% load([ResultPath 'cPDist/cPDistMatrix.mat']);
% % figure;
% % imagesc(cPDistMatrix./max(cPDistMatrix(:))*64);
% % axis equal;
% % axis([1,GroupSize,1,GroupSize]);
% 
% load([ResultPath 'cPDist/cPlmkMSEMatrix.mat']);
% % figure;
% % imagesc(cPlmkMSEMatrix./max(cPlmkMSEMatrix(:))*64);
% % axis equal;
% % axis([1,GroupSize,1,GroupSize]);
% % 
% cPMSTFeatureFixOff = load([ResultPath 'cPMST/FeatureFixOff/cPMSTDistMatrix.mat']);
% cPMSTFeatureFixOff = cPMSTFeatureFixOff.ImprDistMatrix;
% cPMSTFFofflmkMSEMatrix = load([ResultPath 'cPMST/FeatureFixOff/cPMSTlmkMSEMatrix.mat']);
% cPMSTFFofflmkMSEMatrix = cPMSTFFofflmkMSEMatrix.lmkMSEMatrix;
% cPMSTFeatureFixOn = load([ResultPath 'cPMST/FeatureFixOn/cPMSTDistMatrix.mat']);
% cPMSTFeatureFixOn = cPMSTFeatureFixOn.ImprDistMatrix;
% cPMSTFFonlmkMSEMatrix = load([ResultPath 'cPMST/FeatureFixOn/cPMSTlmkMSEMatrix.mat']);
% cPMSTFFonlmkMSEMatrix = cPMSTFFonlmkMSEMatrix.lmkMSEMatrix;
% 
% cPViterbiFeatureFixOff = load('./results/Teeth/cPViterbi/FeatureFixOff/cPViterbiDistMatrix.mat');
% cPViterbiFeatureFixOff = cPViterbiFeatureFixOff.ImprDistMatrix;
% cPViterbiFFofflmkMSEMatrix = load('./results/Teeth/cPViterbi/FeatureFixOff/cPViterbilmkMSEMatrix.mat');
% cPViterbiFFofflmkMSEMatrix = cPViterbiFFofflmkMSEMatrix.lmkMSEMatrix;
% cPViterbiFeatureFixOn = load('./results/Teeth/cPViterbi/FeatureFixOn/cPViterbiDistMatrix.mat');
% cPViterbiFeatureFixOn = cPViterbiFeatureFixOn.ImprDistMatrix;
% cPViterbiFFonlmkMSEMatrix = load('./results/Teeth/cPViterbi/FeatureFixOn/cPViterbilmkMSEMatrix.mat');
% cPViterbiFFonlmkMSEMatrix = cPViterbiFFonlmkMSEMatrix.lmkMSEMatrix;
% 
% cPViterbiAngleFeatureFixOff = load('./results/Teeth/cPViterbiAngle0.5/FeatureFixOff/cPViterbiDistMatrix.mat');
% cPViterbiAngleFeatureFixOff = cPViterbiAngleFeatureFixOff.ImprDistMatrix;
% cPViterbiAngleFFofflmkMSEMatrix = load('./results/Teeth/cPViterbiAngle0.5/FeatureFixOff/cPViterbilmkMSEMatrix.mat');
% cPViterbiAngleFFofflmkMSEMatrix = cPViterbiAngleFFofflmkMSEMatrix.lmkMSEMatrix;
% cPViterbiAngleFeatureFixOn = load('./results/Teeth/cPViterbiAngle0.5/FeatureFixOn/cPViterbiDistMatrix.mat');
% cPViterbiAngleFeatureFixOn = cPViterbiAngleFeatureFixOn.ImprDistMatrix;
% cPViterbiAngleFFonlmkMSEMatrix = load('./results/Teeth/cPViterbiAngle0.5/FeatureFixOn/cPViterbilmkMSEMatrix.mat');
% cPViterbiAngleFFonlmkMSEMatrix = cPViterbiAngleFFonlmkMSEMatrix.lmkMSEMatrix;
% 
% cPLASTFeatureFixOff = load('./results/Teeth/cPLAST/FeatureFixOff/cPLASTDistMatrix.mat');
% cPLASTFeatureFixOff = cPLASTFeatureFixOff.ImprDistMatrix;
% cPLASTFFofflmkMSEMatrix = load('./results/Teeth/cPLAST/FeatureFixOff/cPLASTlmkMSEMatrix.mat');
% cPLASTFFofflmkMSEMatrix = cPLASTFFofflmkMSEMatrix.lmkMSEMatrix;
% cPLASTFeatureFixOn = load('./results/Teeth/cPLAST/FeatureFixOn/cPLASTDistMatrix.mat');
% cPLASTFeatureFixOn = cPLASTFeatureFixOn.ImprDistMatrix;
% cPLASTFFonlmkMSEMatrix = load('./results/Teeth/cPLAST/FeatureFixOn/cPLASTlmkMSEMatrix.mat');
% cPLASTFFonlmkMSEMatrix = cPLASTFFonlmkMSEMatrix.lmkMSEMatrix;
% 
% figure;
% scatter(cPDistMatrix(:),cPMSTFeatureFixOff(:),10,'g');
% axis equal;
% hold on;
% title('cP distances vs cPMST distances');
% plot([0,0.2],[0,0.2],'r');
% axis([0,0.2,0,0.2]);
% 
% figure;
% scatter(cPlmkMSEMatrix(:),cPMSTFFofflmkMSEMatrix(:),10,'b');
% axis equal;
% hold on;
% title('cP landmark MSEs vs cPMST landmark MSEs');
% plot([0,0.7],[0,0.7],'r');
% axis([0,0.7,0,0.7]);
% 
% [~,colIdx] = max(max(cPlmkMSEMatrix-cPMSTFFofflmkMSEMatrix));
% [~,rowIdx] = max(cPlmkMSEMatrix(:,colIdx)-cPMSTFFofflmkMSEMatrix(:,colIdx));

% % close all;
% % figure;
% % scatter(cPDistMatrix(:),ImprDistMatrix(:),10,'g');
% % axis equal;
% % hold on;
% % title('distances before/after MST improvement');
% % plot([0,0.2],[0,0.2],'r');
% % axis([0,0.2,0,0.2]);
% % 
% % figure;
% % scatter(cPlmkMSEMatrix(:),cPMSTlmkMSEMatrix(:),10,'b');
% % axis equal;
% % hold on;
% % title('landmark MSEs before/after MST improvement');
% % plot([0,0.7],[0,0.7],'r');
% % axis([0,0.7,0,0.7]);

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% compare face mesh downsamplings
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SamplePath = [pwd '/samples/Faces/'];
% MeshName = 'quadraticEdgeCollapseDownsampled';
% G = Mesh('off',[SamplePath MeshName '.off']);
% [F,Inds] = unique(sort(G.F',2),'rows','first');
% if ~isempty(setdiff(1:G.nF,Inds))
%     disp('This mesh contains duplicate faces!');
%     V = G.V;
%     F(setdiff(1:G.nF,Inds),:) = [];
%     clear G;
%     G = Mesh('VF',V,F');
% end
% options.display = 'on';
% options.exclude_boundary = 1;
% dVInds = G.DeleteIsolatedVertex(options);
% if ~isempty(dVInds);
%     disp('This mesh contains non-boundary isolated vertex!');
% end
% dVInds = G.FindBoundaries;
% G.DeleteVertex(dVInds);
% G.DeleteIsolatedVertex;
% 
% G.Centralize('ScaleArea');
% % options.SmoothCurvatureFields = 3;
% % options.DensityLocalWidth = 3;
% % options.ConfMaxLocalWidth = 3;
% % options.GaussMaxLocalWidth = 3;
% % options.GaussMinLocalWidth = 3;
% % options.ADMaxLocalWidth = 3;
% % options.ExcludeBoundary = 1;
% % options.Display = 'off';
% G.ComputeMidEdgeUniformization(options);
% 
% save([SamplePath MeshName '.mat'], 'G');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% explore different effects of local max width
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SamplePath = [pwd '/samples/Faces/'];
% MeshName = 'furthestPointDownsampled';
% % MeshName = 'quadraticEdgeCollapseDownsampled';
% load([SamplePath MeshName '.mat']);
% 
% options.SmoothCurvatureFields = 2;
% options.ConfMaxLocalWidth = 4;
% options.GaussMaxLocalWidth = 6;
% options.GaussMinLocalWidth = 6;
% options.ADMaxLocalWidth = 3;
% options.ExcludeBoundary = 1;
% options.Display = 'off';
% 
% G.ExtractFeatures(options);
% 
% figure;G.draw('ConfMax');
% figure;G.draw('GaussMax');
% figure;G.draw('GaussMin');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Compare cP with cPMST and cPMST_FeatureFixing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
% Names = {'B03','w02'};
% SamplePath = [pwd '/samples/Teeth/'];
% cPMaps_path = '/media/trgao10/Work/MATLAB/ArchivedResults/Teeth/cPDist/cPMapsMatrix.mat';
% cPDist_path = '/media/trgao10/Work/MATLAB/ArchivedResults/Teeth/cPDist/cPDistMatrix.mat';
% cPMSTMaps_path = '/media/trgao10/Work/MATLAB/ArchivedResults/Teeth/cPMST/FeatureFixOff/cPMSTMapsMatrix.mat';
% cPMSTFeatureFixMaps_path = '/media/trgao10/Work/MATLAB/ArchivedResults/Teeth/cPMST/FeatureFixOn/cPMSTMapsMatrix.mat';
% data_path = '~/Work/MATLAB/DATA/PNAS/';
% 
% taxa_code = load([data_path 'teeth_taxa_table.mat']);
% taxa_code = taxa_code.taxa_code;
% TAXAind = cellfun(@(name) find(strcmpi(taxa_code,name)),Names);
% 
% disp('loading all cPdist...');
% load(cPDist_path);
% disp('loaded');
% 
% disp('loading all cPmaps...');
% load(cPMaps_path); % load cell array "cPmapsMatrix"
% cPmaps = {cPMapsMatrix{TAXAind(1),TAXAind(2)},cPMapsMatrix{TAXAind(2),TAXAind(1)}};
% disp('loaded');
% 
% disp('loading all cPMSTmaps...');
% load(cPMSTMaps_path); % load cell array "ImprMapsMatrix"
% cPMSTmaps = {ImprMapsMatrix{TAXAind(1),TAXAind(2)},ImprMapsMatrix{TAXAind(2),TAXAind(1)}};
% disp('loaded');
% 
% disp('loading all cPMSTFeatureFixmaps...');
% load(cPMSTFeatureFixMaps_path); % load cell array "ImprMapsMatrix"
% FeatCPMSTmaps = {ImprMapsMatrix{TAXAind(1),TAXAind(2)},ImprMapsMatrix{TAXAind(2),TAXAind(1)}};
% disp('loaded');
% 
% 
% GM = load([SamplePath Names{1} '.mat']);
% GM = GM.G;
% GN = load([SamplePath Names{2} '.mat']);
% GN = GN.G;
% 
% %%% original cP maps
% options.LandmarksPath = [data_path 'landmarks_teeth.mat'];
% options.MeshesPath = [data_path 'meshes/'];
% options.landmarks = 'on';
% options.type = 'full';
% options.ShowCPValue = 'on';
% 
% [~,R,~] = MapToDist(GM.V,GN.V,cPMapsMatrix{TAXAind(1),TAXAind(2)},GM.Aux.VertArea);
% tGM = Mesh(GM);
% tGM.V = R*GM.V;
% 
% ViewTeethCPvsCPMSTvsFeatCPMST(tGM,GN,cPmaps,cPMSTmaps,FeatCPMSTmaps,options);




