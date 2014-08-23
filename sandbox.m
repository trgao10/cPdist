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

GroupSize = 116;

load('./results/Teeth/cPdist/cPdistMatrix.mat');
figure;
imagesc(cPdistMatrix./max(cPdistMatrix(:))*64);
axis equal;
axis([1,GroupSize,1,GroupSize]);

cPlmkMSEMatrix = load('./results/Teeth/cPdist/lmkMSEMatrix.mat');
cPlmkMSEMatrix = cPlmkMSEMatrix.lmkMSEMatrix;
figure;
imagesc(cPlmkMSEMatrix./max(cPlmkMSEMatrix(:))*64);
axis equal;
axis([1,GroupSize,1,GroupSize]);

load('./results/Teeth/cPMST/FeatureFixOff/cPMSTDistMatrix.mat');
figure;
imagesc(ImprDistMatrix./max(cPdistMatrix(:))*64);
axis equal;
axis([1,GroupSize,1,GroupSize]);

cPMSTlmkMSEMatrix = load('./results/Teeth/cPMST/FeatureFixOff/cPMSTlmkMSEMatrix.mat');
cPMSTlmkMSEMatrix = cPMSTlmkMSEMatrix.lmkMSEMatrix;
figure;
imagesc(cPMSTlmkMSEMatrix./max(cPMSTlmkMSEMatrix(:))*64);
axis equal;
axis([1,GroupSize,1,GroupSize]);

close all;

figure;
scatter(cPdistMatrix(:),ImprDistMatrix(:),10,'g');
axis equal;
hold on;
title('distances before/after MST improvement');
plot([0,0.2],[0,0.2],'r');
axis([0,0.2,0,0.2]);

figure;
scatter(cPlmkMSEMatrix(:),cPMSTlmkMSEMatrix(:),10,'b');
axis equal;
hold on;
title('landmark MSEs before/after MST improvement');
plot([0,0.7],[0,0.7],'r');
axis([0,0.7,0,0.7]);




