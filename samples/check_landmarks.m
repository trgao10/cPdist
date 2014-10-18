Name = '04';

A = Mesh('off','~/Work/MATLAB/DATA/Clement/meshes/01.off');
load('~/Work/MATLAB/DATA/Clement/landmarks_clement.mat');
landmarks = squeeze(PP(4,:,:));

A.draw();
hold on;
scatter3(landmarks(:,1),landmarks(:,2),landmarks(:,3),30,'b','filled');

