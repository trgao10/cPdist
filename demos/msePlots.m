load('/media/trgao10/Work/MATLAB/ArchivedResults/PNAS/cPDist/cPlmkMSEMatrix.mat');
load('/media/trgao10/Work/MATLAB/ArchivedResults/PNAS/cPViterbiAngle0.25/FeatureFixOff/cPViterbilmkMSEMatrix.mat')

plot(cPlmkMSEMatrix(:),lmkMSEMatrix(:),'b.');
axis equal;
ylim([0,0.75]);
xlim([0,0.75]);
hold on;
plot([0,0.75],[0,0.75],'r-');

xlabel('cPDist');
ylabel('cPViterbi-FFoff');